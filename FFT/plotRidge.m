function [] = plotRidge(R,Ev,Ts,N,sampleName,binSize,zeroPadding,derivativeIndicator)
%will plot the maximum-energy time-period ridge extracted from the STFT spectrogram or CWT scalogram
%R: ridge data provided as cell array, with each cell including the data for one sample;
%R{1,i} - time [min]; R{2,i} - periods with maximum energy at a given time point [min]
%Ev:  timings of cellular events of interest in [min] starting from
%nucleoid segregation; provided as a table with labels of event type
%N: number of time points for each sample

outputPath = ['/Users/muxika/MATLAB/result/' sampleName];
mkdir(outputPath);
outputPathSpectr = [outputPath '/ridge_vs_genEv.pdf'];


Nsamples = size(R,2);
Pmax = max(cellfun(@(c) max(c(:,2)), R));
Evnames = Ev.Properties.VariableNames;
NEv = size(Evnames,2); 

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
            counter = counter +1;
            subplot(3,2,counter)
            %plot ridge
            plot(R{j}(:,1),R{j}(:,2))
            ylabel('Period [min]')
            ylim([2*Ts/60 Pmax])
            xlim([0 (Ts/60)*N(j)])
            for l = 2:NEv
                vline(minutes(Ev.(l)(j)),':',Evnames{l},[0 1],{'Rotation', 90},[],[],[2*Ts/60 Pmax/2]);
            end
            xlabel('Nucleoid-cycle time [min]')
            title(['Sample' num2str(j)])
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

end

