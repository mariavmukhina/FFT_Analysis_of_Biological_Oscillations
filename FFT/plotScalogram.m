function [] = plotScalogram(L,Ts,N,DCLimit,sampleName,binSize,zeroPadding,derivativeIndicator)
%will plot continuous wavelet transform (CWT) scalograms
%L : samples in time domain, after detrending or derivative and zero-padding)
%Ts: sampling period [sec]
%N: number of time points in samples before zero padding
%DCLimit: cutoff for long periodicities

outputPath = ['/Users/muxika/MATLAB/result/' sampleName];
mkdir(outputPath);
outputPathSpectr = [outputPath '/scalogram.pdf'];
Nsamples = size(L,2);

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
parameters(6,2) = {num2str(Ts)};

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
            t = ((0:1:N(j)-1).*(Ts/60)); % calculate time axis in min
            [wt,period] = cwt(L{j},'amor',minutes(Ts/60)); % continuous wavelent transform using Morlet wavelet
            tridge = tfridge(abs(wt(:,1:N(j))),minutes(period)); % find the ridge (in min)
            counter = counter +1;
            subplot(3,2,counter)
            %plot scalogram
            surface(t,minutes(period),abs(wt(:,1:N(j))))
            axis tight
            shading flat
            ylabel('Period [min]')
            ylim([2*Ts/60 DCLimit])
            xlabel('Nucleoid-cycle time [min]')
            title(['Sample' num2str(j)])
            ridge{j} = [t tridge']; % save ridge (1) time in min (2) period in min
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

save([outputPath '/CWTridge.mat'],'ridge');
end

