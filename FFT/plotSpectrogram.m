function [] = plotSpectrogram(L,Fs,nfft,N,Pmax,DCLimit,sampleName,binSize,zeroPadding,derivativeIndicator)
%will plot short time Fouirier transform (STFT) spectrograms
%L : samples in time domain, after detrending or derivative and zero-padding
%Fs: sampling rate [mHz]
%nfft: number of FFT points (length of padded sample)
%N: number of samples (before zero-padding)
%Pmax: vector if the window is set for each sample (e.g. period of oscillation with maximum intensity) or a single value [min] 

outputPath = ['/Users/muxika/MATLAB/result/' sampleName];
mkdir(outputPath);
outputPathSpectr = [outputPath '/spectrogram.pdf'];
Nsamples = size(L,2);
if numel(Pmax) == 1
    Pmax = Pmax.*ones(Nsamples,1); % [min] for constant window
end
window = round(Pmax./(1000/Fs/60)); % convert periods in min to number of points

%create report object
import mlreportgen.dom.*
import mlreportgen.report.*
rpt = Report(outputPathSpectr,'pdf'); 

%create table with FFT parameters for the whole dataset
parameters = cell(6,2);
parameters(1,1) = {'FFT parameters'};
parameters(2:6,1) = {'dataset'; 'binSize, in mHz'; 'zero padding'; 'FFT on derivative'; 'sampling period, sec   '};
parameters(2,2) = {sampleName};
if size(binSize,2) ~= 1
    minBin = min(binSize);
    maxBin = max(binSize);
    binRange = {['bin size is in the range between ' num2str(minBin) ' and ' num2str(maxBin)]};
    parameters(3,2) = binRange;
else
    parameters(3,2) = {num2str(binSize)};
end
parameters(4,2) = {num2str(zeroPadding)};
parameters(5,2) = {num2str(derivativeIndicator)};
parameters(6,2) = {num2str(1000/Fs)};

table = Table(parameters);
add(rpt,table);

for k = 1:6:Nsamples
    n1              = k;
    n2              = k+5;
    counter         = 0;
    fig = figure('visible','off');
    hold on;
    try
        for j = n1:n2
            [~,f,t,psd] = spectrogram(L{j},window(j),window(j)-1,nfft,Fs*1e-3,'yaxis'); % calculate spectrogram in Hz
            t = t./60; % convert time vector to min
            f = f.*1000; %convert frequency vector to mHz
            DC = numel(find(f<=DCLimit)); % remove slow frequencies
            f(1:DC) = [];
            psd(1:DC,:) = [];
            fridge = tfridge(psd,f); % find the ridge
            [X,Y] = meshgrid(t,f);
            counter = counter +1;
            subplot(3,2,counter)
            %plot spectrogram
            surf(X,Y,psd)
            view(2)
            shading interp
            axis tight
            ylabel('Frequency [mHz]')
            xlabel('Nucleoid-cycle time [min]')
            xlim([0+window(j) N(j)*(1000/60/Fs)])
            title(['Sample' num2str(j) '; STFT window, min ' num2str(window(j).*(1000/Fs/60))])
            z = find(fridge); % remove zero elements due to zero padding
            ridge{j} = [t(z)'-t(z(1)) (1000/60)./fridge(z)]; % save ridge (1) time in min (2) period in min
        end
    catch
        %to prevent stopping because of the error occuring when number of
        %samples is not multiple of 6
    end
figReporter = Figure(fig);
figImg = Image(getSnapshotImage(figReporter,rpt));
figImg.Style = [figImg.Style {ScaleToFit}];
add(rpt,figImg)

close gcf             
end
    
close(rpt);

save([outputPath '/STFTridge.mat'],'ridge');
end

