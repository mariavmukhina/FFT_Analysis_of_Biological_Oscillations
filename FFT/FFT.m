function [] = FFT(sampPeriod,binSize,zeroPadd,derivative,outputs)
%FFT loads data from .mat files in dataFolder, calculates FFT and generates
%outputs based on input parameters
%
% input .mat files in dataFolder; can be either num arrays arranged so column represents a
%            single nucleoid length curve or cell arrays a(n,1) where each row
%            contains a single nucleoid length curve
% binSize: in mHz, the smallest desirable frequency difference; set if use
% zero padding
%
% sampPeriod: temporal resolution of image acquisition in sec, can be a
%             vector if datasets are taken with different temporal
%             resolution [sec]
%
% zeroPadd: 1 - zero padding, 2 - no zero padding, does not improve the
%           resolution and produces side lobes but provide more accurate amplitude estimation
% derivative: 1 - use signal derivative to calculate FFT ("velocity of
%           displacement", small fast oscillations have higher energy), 0 - use
%           signal to calculate FFT ("amplitude of
%           displacement", large slow oscillations have higher energy)
% outputs: [1 1 1 1 0.6 1 1 1 1 1 1 3] - vector of zeros and ones, if (outputs(i) == 1, generates: 
%                      1) induvidual PSD spectra; 
%                      2) avg PSD for all datasets (select only for zero-padded spectra);
%                      3) 2D energy vs period histogram;
%                      4) CDF histograms for all datasets (energy and period of the most intense FFT peaks);
%                      5) kernel smoothing function estimates for energy and periods of the most intense FFT peaks, NB! outputs(x) can be ~=1, where x
%                      defines correction factor for kernel bandwidth, <1 decreases smoothing, >1 increases smoothing;
%                      6) save all intense FFT frequencies (as defined by findpeak function) and corresponding periods and powers to .mat files;
%                      7) plot samples statistics (number of samples in dataset and distributions of sample lengths before zero-padding);
%                      8) 3D PSD surfaces representing PSD spectra for each sample in datasets;
%                      9) induvidual length curves filtered with low-pass FFT filter;
%                     10) time-frequency spectrograms (STFT) for each sample in datasets;
%                     11) Continuous Wavelet Transform (CWT) for each sample in datasets (requires wavelet toolbox).
%                     12) plot difference spectra for pairs(n0,n) of datasets where n0 is the dataset which is being subtructed, typically WT slow growth; put no to outputs vector
%  
% By Maria Mukhina for Kleckner Lab.
% mmukhina@fas.harvard.edu
% 24.02.2020
% GNU GENERAL PUBLIC LICENSE v3


dataFolder                = '/Users/muxika/MATLAB/temp';
outputFolder              = '/Users/muxika/MATLAB/result';
dataPath                  = fullfile(dataFolder, '*.mat');

% structure with files list
files = dir(dataPath); 
[~, idx] = natsortfiles({files.name});
files = files(idx);
% number of files
numFiles = numel({files.name});
%set colors
C = linspecer(numFiles,'sequential');

F_Nyquist = 1000/(4*max(sampPeriod))-0.03; % in mHz defined by Nyquist theorem, 0.03 is empirically defined to be able to compare different rates
DCLimit = 0.28; % in mHz 0.6 = 27.8 min (for all datasets, 0.28 mHz = 60 min for slow growth) chosen empirically to remove DC components, many samples do not have enough datapoints to detect such long period  

%% calculate sampling rate for each dataset
if length(sampPeriod) ~= 1
    Fs = 1000./sampPeriod; % sampling rate in mHz  
elseif length(sampPeriod) == 1
    Fs(1:numFiles) = 1000/sampPeriod; % sampling rate in mHz  
    sampPeriod(1:numFiles) = sampPeriod;
end

%% cycle through datasets
for n = 1 : numFiles

    
    if zeroPadd == 1  
        % find length of signal necessary to achieve desired bin size
        nfft = round(Fs(n)/binSize);
    end
    
    %% import data
    fileName              = files(n).name; 
    signal                = importdata([files(n).folder '/' fileName]);
    if iscell(signal) == 1
       signal = cellListToArray(signal,2);
    end
    numberOfSamples       = size(signal,2);
    
    %% calculate FFT, PSD, find all FFT peaks 
    for j = 1:numberOfSamples
        oneSample         = signal(:,j);
        oneSample(isnan(oneSample))=0; % substitute to zeros if sample is padded with NANs
        oneSample(oneSample == 0)=[]; % delete all padded zeros
        signalLength = size(oneSample,1);
        avNucLength(j) = mean(oneSample);
        % save lengths of individuals signals before zero padding
        N(j) = signalLength;
        if derivative == 1  
            %differentiate raw signal
            oneSample      = diff(oneSample);
        elseif isempty(derivative) || derivative == 0
            %detrend
            oneSample = subtructBL(oneSample,5,sampPeriod);
        end
        
        %save data after pre-processing
        Lcurves(j) = {oneSample};
        
        if zeroPadd == 1  
            pad = nfft - signalLength +1;
            if pad < 0
                warning('Bin size is too small')
            end
            % zero padding to paddedSignalLength
            oneSample = [oneSample;zeros(pad,1)];
            % save padded length curves
            Lcurves_pad(j) = {oneSample};
            signalLength = size(oneSample,1); %FFT length for desired resolution in frequency space
            % find frequencies to calculate Fourier for zero-padded signal
            f  = binSize*(0:signalLength/2-1); % in mHz

        elseif isempty(zeroPadd) || zeroPadd == 0
            binSize(j) = Fs(n)/signalLength;
            % find frequencies to calculate Fourier for non-padded signal
            f  = binSize(j)*(0:signalLength/2-1); % in mHz
        end  
        %normalize amplitudes of the length curves to [0 1] ******** optional %************
        oneSample = oneSample./max(abs(oneSample));
        %caculate Fourier spectrum of the signal
        FFT          = fft(oneSample,signalLength); 
        %calculate equivalent noise bandwidth assuming there is no windowing (or rectangular window is applied which is equivalent)
        ENBW = (Fs(n)/min(Fs))/N(j); % 
        % measure energy of various frequencies (power spectral dencity), *2 to get single-sided spectrum
        PSD(j)        = {(2*abs(FFT)*ENBW).^2}; 
        %Select Freq within the Nyquist interval [mHz]
        NyquistInterval            = (DCLimit <= f & f <= F_Nyquist);
        % remove all freq outside of Nyquist interval
        f(~any(NyquistInterval,1)) = [];
        f_Nqst(j) = {f};
        P_Nqst(j)  = {PSD{j}(1:floor(signalLength/2),1)};
        P_Nqst{j}(~any(NyquistInterval,1)) = [];
        %find all frequencies and their energy
        [Energies{j},locs{j}] = findpeaks(P_Nqst{j}); % 'MinPeakProminence', 0.5
        Frequencies{j} = f_Nqst{j}(locs{j});
        Periods{j} = 16.6666./Frequencies{j}; %min
          
    end

%% Prepare data for optional outputs (calculate averages, find max, etc) %%%%%%%%%%%%%%%%%%
%% Prepare data for output 2
if (outputs(2) == 1 || outputs(12) ~=0) & zeroPadd == 1
    % Save Fourier frequencies and PSD amplitude for future use
    FourierFreq_allDataSets(n)  = {f_Nqst};
    PSD_allDataSets(n)  = {P_Nqst};
    %Calculate averaged Power Spectral Density spectrum
    PSDaverage = mean(cell2mat(P_Nqst),2);
    PSDavg_allDataSets(n) = {PSDaverage/j}; 
elseif (outputs(2) == 1 || outputs(12) ~=0) & (zeroPadd == 0 || isempty(zeroPadd))
    disp('Need to use zero padding to calculate average PSD spectrum')
end

%% Save .mat outputs for Energies,Frequencies,Periods %%%%%%%%%%%%%%%%%%%%%
fileName = extractBefore(fileName,'.mat');
outputFileName(n) = {strrep(fileName,'_',' ')};

if outputs(1,6) == 1
    PSDoutput.Energies = Energies;
    PSDoutput.Frequencies = Frequencies;
    PSDoutput.Periods = Periods;
    
    if isempty(derivative) || derivative == 0
        save([outputFolder '/' char(fileName) '_PSDpeaks.mat'],'PSDoutput');
    elseif derivative == 1
        save([outputFolder '/' char(fileName) '_PSDonDerivPeaks.mat'],'PSDoutput');
    end
end
%% Prepare data for output 4 & 5 & 7 & 8 & 10 & 12
if outputs(4) == 1 || outputs(5) ~= 0 || outputs(8) == 1 || outputs(7) == 1 || outputs(10) == 1 || outputs(12) ~= 0
    %Find FFT peak with maximum intensity
    EnergyMAX = cellfun(@(c) max(c), Energies);
    I_period = cellfun(@(c) find(c == max(c)), Energies);
    EnergyMAX_allDataSets(n) = {EnergyMAX};
    EnergyMAX_inDataset(n) = max(EnergyMAX);
    for t = 1:size(Periods,2)
        PeriodMAX(t) = Periods{1,t}(I_period(t));
    end
    PeriodMAX_allDataSets(n) = {PeriodMAX};
    PeriodMAX_inDataset(n) = max(PeriodMAX);   
end


%% Prepare data for sample statistics (output 7)
if outputs(7) == 1
    %save average nucleoid length for each cell in dataset
    avNucL_allDataSets(n) = {avNucLength};
    actualSampLength_allDataSets(n) = {N};
    numSamp_allDataSets(n) = {num2str(j)};
end

%% Plot Individual PSD spectra (output 1) %%%%%%%%%%%%%%%%%%%%%%%%%

if outputs(1,1) == 1
    % Plot induvidual PSD for all datasets
    plotIndividualPSDspectra(P_Nqst,f_Nqst,Lcurves,N*sampPeriod(n)/60,fileName,binSize,zeroPadd,derivative,sampPeriod(n),DCLimit);
end

%% Plot 2D Energy Histogram (output 3) %%%%%%%%%%%%%%%%%%%%%%%%%    
if outputs(1,3) == 1
    %Plot 2D energy histogram
    plot2DEnergyHistogram(Energies,Periods);
end

%% Plot 3D PSD surface (output 8)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if outputs(8) == 1 & zeroPadd == 1
    plotPSDsurf(P_Nqst,f_Nqst{1},'P',PeriodMAX,fileName,sampPeriod(n),binSize,zeroPadd,derivative); % 'L' - Lcurves; 'P' - PeriodMAX
elseif outputs(8) == 1 & zeroPadd == 0
    disp('you need to zero-padd in order to plot 3D PSD surface for the whole dataset')
end

%% Plot filtered length curves (output 9) %%%%%%%%%%%%%%%%%%%%%%%%%

if outputs(1,9) == 1 & zeroPadd == 1
    % Filter high frequencies out and plot induvidual length curves for all datasets 
	cutoff = 5; %min
    plotFilteredLengthCurves(Lcurves_pad,cutoff,Fs(n),nfft,N,fileName,binSize,zeroPadd,derivative);
elseif outputs(1,9) == 1 & zeroPadd == 0
    disp('you need to zero-padd in order to plot length curves filtered with low-pass filter')
end

%% Plot STFT spectrograms (output 10) %%%%%%%%%%%%%%%%%%%%%%%%%

if outputs(1,10) == 1 & zeroPadd == 1
    % Calculatetes and plots short time fourier transform (STFT) spectrograms and their max energy ringes for induvidual length curves for all datasets 
    plotSpectrogram(Lcurves_pad,Fs(n),nfft,N,10,DCLimit,fileName,binSize,zeroPadd,derivative); % 
elseif outputs(1,10) == 1 & zeroPadd == 0
    disp('you need to zero-padd in order to plot spectrograms')
end

%% Plot CWT scalograms (output 11) %%%%%%%%%%%%%%%%%%%%%%%%%

if outputs(1,11) == 1 & zeroPadd == 1
    % % Calculatetes and plots continuous wavelet transform (CWT) scalograms and their max energy ringes for induvidual length curves for all datasets 
    plotScalogram(Lcurves_pad,sampPeriod(n),N,(1000/60)./DCLimit,fileName,binSize,zeroPadd,derivative);
elseif outputs(1,11) == 1 & zeroPadd == 0
    disp('you need to zero-padd in order to plot spectrograms')
end

%% clear all variables in j cycle
clear N Lcurves p f_Nqst P_Nqst Energies Frequencies Periods locs ENBW EnergyMAX PeriodMAX I X Y

if size(binSize,2) ~= 1
    binSize = [];
end

end

%% Finish preparation of the data for output 4 & 5
if outputs(4) == 1 || outputs(5) ~= 0 
    % pick dataset with highest values for power/period
    
    m = find(EnergyMAX_inDataset == max(EnergyMAX_inDataset));
    mm = find(PeriodMAX_inDataset == max(PeriodMAX_inDataset));
    
    %find bin edges for Emax and Pmax histograms and rug plots
    [~,edges_E] = histcounts(EnergyMAX_allDataSets{m(1)},50);
    [~,edges_P] = histcounts(PeriodMAX_allDataSets{mm(1)},50);

    %find x points for kernel smoothing function estimates for energy and periods distributions
    [~,xi_E,bw_E] = ksdensity(EnergyMAX_allDataSets{m(1)},'Kernel','epanechnikov'); 
    [~,xi_P,bw_P] = ksdensity(PeriodMAX_allDataSets{mm(1)},'Kernel','epanechnikov');

    % correct kernel bandwidth, choose < 1 if kernel smoothing function is too wide
    bw_P = bw_P*outputs(5);
    bw_E = bw_E*outputs(5);

    %calculate adder for jitter in rug plot
    adder_E = (edges_E(2) - edges_E(1))/4/numFiles;
    adder_P = (edges_P(2) - edges_P(1))/4/numFiles;

    for i = 1:n
        %Calculate kernel smoothing function estimates for energy and periods distributions (output 5)
        [PDF_E_allDataSets{i},~] = ksdensity(EnergyMAX_allDataSets{i},xi_E,'Kernel','epanechnikov','Bandwidth',bw_E); 
        [PDF_P_allDataSets{i},~] = ksdensity(PeriodMAX_allDataSets{i},xi_P,'Kernel','epanechnikov','Bandwidth',bw_P); 

        %Calculate histogram bin counts for rug plot (output 5)
        [binCountsPDF_E_allDataSets{i},~] = histcounts(EnergyMAX_allDataSets{i},edges_E,'Normalization','pdf');
        [binCountsPDF_P_allDataSets{i},~] = histcounts(PeriodMAX_allDataSets{i},edges_P,'Normalization','pdf');

        %Calculate histogram bin counts for CDF linear histograms (output 4)
        [binCountsCDF_E_allDataSets{i},edgesCDF_E_allDataSets{i}] = histcounts(EnergyMAX_allDataSets{i},50,'Normalization','cdf');
        [binCountsCDF_P_allDataSets{i},edgesCDF_P_allDataSets{i}] = histcounts(PeriodMAX_allDataSets{i},50,'Normalization','cdf');
        
    end
end

%% Plot CDF histograms for max energies and periods for all datasets %%%%%%%%

if outputs(1,4) == 1
    
    %%%%%%%%%%%Plot energy histogram
    figure; 
    ax1 = axes;
    set(ax1,'NextPlot','replacechildren', 'ColorOrder',C,'Box','on');
    ax1.YLabel.String ='CDF';
    ax1.XLabel.String ='Power of the most intense peak';
    if derivative == 1
        ax1.XLim = [0 6];
    elseif derivative == 0
        ax1.XLim = [0 6];
    end

    for i = 1:n
        %hold(ax1,'on')
        plot(ax1,edgesCDF_E_allDataSets{i}(1,1:end-1),binCountsCDF_E_allDataSets{i},'LineWidth', 2);
        hold on;
        %hold(ax1,'off')
    end
    legend(ax1,outputFileName,'Location','east');
    title({['binSize ', num2str(binSize), ' mHz; zero Padding ', num2str(zeroPadd), '; FFT on derivative ', num2str(derivative), ';'];[ 'sampPeriod ', num2str(sampPeriod),' sec']})
    
    %%%%%%%%%%%Plot period histogram
    figure; 

    ax2 = axes('NextPlot','replacechildren', 'ColorOrder',C,'Box','on'); 
    ax2.YLabel.String ='CDF';
    ax2.XLabel.String ='Period of the most intense peak [min]';
    
    if derivative == 1
        ax2.XLim = [2 12];
    elseif derivative == 0
        ax2.XLim = [2 27];
    end

    for i = 1:n
        plot(ax2,edgesCDF_P_allDataSets{i}(1,1:end-1),binCountsCDF_P_allDataSets{i},'LineWidth', 2);
        hold on;
    end
    legend(ax2,outputFileName,'Location','east');
    title({['binSize ', num2str(binSize), ' mHz; zero Padding ', num2str(zeroPadd), '; FFT on derivative ', num2str(derivative), ';'];[ 'sampPeriod ', num2str(sampPeriod),' sec']})
end

%% Plot kernel smoothing function estimates for max energy and periods distributions for all datasets
if outputs(5) ~= 0 
    %energy
    figure; 
    axes('NextPlot','replacechildren', 'ColorOrder',repelem(C,2,1));
    jitter_E = 1; % jitter is necessary to see all the distributions
    for i = 1:n
        plot(xi_E,PDF_E_allDataSets{i},'LineWidth', 2,'DisplayName',outputFileName{i});
        hold on
        stem(edges_E(1,1:end-1)*jitter_E,0.2*binCountsPDF_E_allDataSets{i},'Marker','none','HandleVisibility','off');
        hold on
        jitter_E = jitter_E + adder_E; 
    end
    if derivative == 1
        xlim([0 8])
    elseif derivative == 0
        xlim([0 5])
    end
    ylabel('PDE')
    xlabel('Power of the most intense peak')
    legend('show');
    title({['binSize ', num2str(binSize), ' mHz; zero Padding ', num2str(zeroPadd), '; FFT on derivative ', num2str(derivative), ';'];[ 'sampPeriod ', num2str(sampPeriod),' sec']})

    %period
    figure; 
    axes('NextPlot','replacechildren', 'ColorOrder',repelem(C,2,1));
    jitter_P = 1;
    for i = 1:n
        plot(xi_P,PDF_P_allDataSets{i},'LineWidth', 2,'DisplayName',outputFileName{i});
        hold on
        stem(edges_P(1,1:end-1)*jitter_P,0.2*binCountsPDF_P_allDataSets{i},'Marker','none','HandleVisibility','off');
        hold on
        jitter_P = jitter_P + adder_P;
    end
    if derivative == 1
        xlim([2 7])
    elseif derivative == 0
        xlim([2 27])
    end
    ylabel('PDE')
    xlabel('Period of the most intense peak [min]')
    legend('show');
    title({['binSize ', num2str(binSize), ' mHz; zero Padding ', num2str(zeroPadd), '; FFT on derivative ', num2str(derivative), ';'];[ 'sampPeriod ', num2str(sampPeriod),' sec']})

end


%% Plot average PSD for all datasets
if outputs(2) == 1
    
    figure; 
    axes('NextPlot','replacechildren', 'ColorOrder',C);
    for i = 1:n
        plot(FourierFreq_allDataSets{i}{1},PSDavg_allDataSets{i}); % take first cell because frequencies are the same for zero-padded spectra
        hold on;
    end
    xlabel('Frequency [mHz]')
    ylabel({'Average Power Spectral Density'; 'normalized by number of samples in dataset[(px/s)^2*mHz]'})
    legend(outputFileName);
end

%% Plot samples statistics for all datasets
if outputs(7) == 1
    %calculate distribution of average volumes of nucleoids assuming that nucleoid is a perfect cilinder
    for i = 1:n
        nucVol{i} = pi*(avNucL_allDataSets{i}./5).^2.*avNucL_allDataSets{i};
        [binCountsCDF_avNucL_allDataSets{i},edgesCDF_avNucL_allDataSets{i}] = histcounts(nucVol{i},50,'Normalization','cdf'); 
    end
    
    %plot distribution of average lengths of nucleoids
    figure; 
    axes('NextPlot','replacechildren', 'ColorOrder',C);
    for i = 1:n
        plot(edgesCDF_avNucL_allDataSets{i}(1,1:end-1),binCountsCDF_avNucL_allDataSets{i},'LineWidth',2); 
        hold on;
    end
    xlabel('Average nucleoid length [px]')
    ylabel('Number of samples')
    legend(outputFileName); 
end

%% Plot difference PSD spectra
if outputs(12) ~=0
    plotDifferencePSD(PSD_allDataSets,FourierFreq_allDataSets,PeriodMAX_allDataSets,outputs(12),{files.name},sampPeriod,binSize,zeroPadd,derivative);
end

