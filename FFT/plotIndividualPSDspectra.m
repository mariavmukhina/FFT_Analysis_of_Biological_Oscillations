function [] = plotIndividualPSDspectra(PSD,f,lengthCurves,sampleLengths,sampleName,binSize,zeroPadding,derivativeIndicator,sampPeriod,DClimit)

%power: FFT power within Nyquist interval
%f:     Fourier frequencies within this interval
%lengthCurves: samples in time domain, displacement or derivative before zero padding
%sampleLengths: lengths of the samples before zero-padding in min

tmin = (1000/60)/8.33;
tmax = (1000/60)/DClimit;

outputPath = ['/Users/muxika/MATLAB/result/' sampleName];
mkdir(outputPath);
outputPath = [outputPath '/PSDspectra.pdf'];

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
parameters(6,2) = {num2str(sampPeriod)};

table = Table(parameters);
add(rpt,table);

%plot figures of individual PSD spectra and corresponding length curves in
%time domain

Nsamples = size(PSD,2);

PSD_max = max(cellfun(@(c) max(c), PSD));
PSD_max = PSD_max*0.7; % to set the y limit at 80% of the max value

length_max = max(cellfun(@(c) max(c), lengthCurves));
length_N = max(cellfun(@(c) size(c,1), lengthCurves));

for k = 1:3:Nsamples
    n1              = k;
    n2              = k+2;
    counter         = 0;
    fig = figure('visible','off');
    hold on;
    try
        for j = n1:n2
            %plot PSD vs f
            counter = counter +1;
            subplot(3,2,counter)
            plot((1000/60)./f{j},PSD{j}) % covert frequencies in mHz to min 
            %xlabel('Frequency [mHz]')
            xlabel('Period [min]')
            ylabel('Normalized PSD')
            %ylim([0 PSD_max]);
            %xlim([0.6 8.33]);
            xlim([tmin tmax]);
            title(['Sample' num2str(j) ]) 
            counter = counter +1;
            %plot length curves
            subplot(3,2,counter)
            plot(lengthCurves{j})
            xlabel('Frame')
            ylabel('Length [px]')
            xlim([0 length_N])
            ylim([(-1)*length_max length_max])
            if derivativeIndicator == 0
                title(['Length curve, acquisition duration ' num2str(sampleLengths(j)) ' min'])
            elseif derivativeIndicator == 1
                title(['Length curve derivative, acquisition duration ' num2str(sampleLengths(j)) ' min'])
            end
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

