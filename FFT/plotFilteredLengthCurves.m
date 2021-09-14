function [] = plotFilteredLengthCurves(L,cutoff,Fs,nfft,N,sampleName,binSize,zeroPadding,derivativeIndicator)

%L : samples in time domain, after detrending or derivative and after zero padding
%Fs: sampling rate [mHz]
%nfft: number of FFT points (legth of padded sample)
%N: number of samples before zero-padding


outputPath = ['/Users/muxika/MATLAB/result/' sampleName];
mkdir(outputPath);
outputPath = [outputPath '/LengthCurvesLowPassFilter.pdf'];
Nsamples = size(L,2);

%create report object
import mlreportgen.dom.*
import mlreportgen.report.*
rpt = Report(outputPath,'pdf'); 

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

% design low pass filter
Fstop  = 1/(60*cutoff); %Hz cut-off frequency; difference between Fpass and Fstop defines filter transition width
Fpass = Fstop - Fstop/10; %Hz

d = designfilt('lowpassfir','SampleRate',Fs*1e-3,'PassbandFrequency',Fpass, ...
         'StopbandFrequency',Fstop,...
         'StopbandAttenuation',30);
D = round(mean(grpdelay(d))); %find the group delay of the filter 

for i = 1:Nsamples
    % perform filtering
    L_lp{i} = fftfilt(d,[L{i}; zeros(D,1)],nfft);  % append D zeros at the end of the input data to ensure that all data are flushed out of the filter
    % Shift data to compensate for delay
    L_lp{i} = L_lp{i}(D+1:end);                  %
    %remove padding zeros
    L_lp{i}(N(i)+1:end) = [];
    L{i}(N(i)+1:end) = []; 
    %calculate time axis
    t{i}     = ((0:1:N(i)-1).*(1000/Fs/60))';
end


%plot individual length curves in time domain, filtered with low pass
%filter and original


length_max = max(cellfun(@(c) max(c), L));
length_N = max(cellfun(@(c) size(c,1), L));

for k = 1:3:Nsamples
    n1              = k;
    n2              = k+2;
    counter         = 0;
    fig = figure('visible','off');
    hold on;
    try
        for j = n1:n2
            %plot filtered length curves
            counter = counter +1;
            subplot(3,2,counter)
            plot(t{j},L_lp{j}) 
            xlabel('Nucleoid-cycle time [min]')
            ylabel('Length [px]')
            xlim([0 length_N])
            ylim([(-1)*length_max length_max])
            title(['Sample' num2str(j) '; cutoff ' num2str(1/60/Fstop) ' min']) 
            counter = counter +1;
            %plot original length curves
            subplot(3,2,counter)
            plot(t{j},L{j})
            xlabel('Nucleoid-cycle time [min]')
            ylabel('Length [px]')
            xlim([0 length_N])
            ylim([(-1)*length_max length_max])
            title('Unfiltered length curve')
        end
    catch
        %to prevent stopping because of the error occuring when number of
        %samples is not multiple of 3
    end
figReporter = Figure(fig);
figImg = Image(getSnapshotImage(figReporter,rpt));
figImg.Style = [figImg.Style {ScaleToFit}];
add(rpt,figImg)

close gcf             
end
    
close(rpt);

end
