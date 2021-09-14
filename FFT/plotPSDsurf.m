function [] = plotPSDsurf(PSD,f,srt,srtData,datasetName,sampPeriod,binSize,zeroPadding,derivativeIndicator)
%will plot PSD(1/f) spectra for each sample in the dataset as 3D surface
%where x axis is a period (1/f), y axis is the numbers of samples in
%dataset, and z(color) axis is PSD
% PSD: PSD amplitudes for each sample in the dataset
% f: Fourier frequencies for each sample in the dataset
% srt: type of soorting to be used to build a PSD surface: 'P' - by period [min] of FFT peak with maximum PSD; 'L' - by length of the nucleoid cycle
% srtData: data based on which is performed sorting in the dataset, can be period [min] of FFT peak with maximum PSD for each sample in the dataset (for 'P' sorting)
% or length of the nucleoid cycles (for 'L' sorting)

outputPath = ['/Users/muxika/MATLAB/result/' datasetName];
mkdir(outputPath);
outputPathSpectr = [outputPath '/datasetPSD.pdf'];

%create report object
import mlreportgen.dom.*
import mlreportgen.report.*
rpt = Report(outputPathSpectr,'pdf'); 

%create table with FFT parameters for the whole dataset
parameters = cell(2,5);
parameters(1,1:5) = {'dataset' 'binSize[mHz]' 'zero_padding' 'FFT_on_derivative' 'sampling_period[sec]'}; 
parameters(2,1) = {datasetName};
parameters(2,2) = {num2str(binSize)};
parameters(2,3) = {num2str(zeroPadding)};
parameters(2,4) = {num2str(derivativeIndicator)};
parameters(2,5) = {num2str(sampPeriod)};
if srt == 'L'
    parameters(3,1) = {'samples sorted by the length of the nucleoid cycle'};
    srtData = cellfun(@(x) size(x,1), srtData);
elseif srt == 'P'
    parameters(3,1) = {'samples sorted by the period of FFT peak with maximum PSD'};
end

table = Table(parameters);
add(rpt,table);

%find number of samples in dataset
s = size(PSD,2);
% greate a grid for plotting
[X,Y] = meshgrid((1000/60)./f,1:s);

%sort samples in dataset by data provided for sorting

[~,idx] = sort(srtData);
PSD = cell2mat(PSD(idx))';

fig = figure('visible','off');
s = surf(X,Y,PSD);
s.EdgeColor = 'none';
c = colorbar;
c.Label.String = 'PSD';
view(-18.698935578363781,72.538181818181812)
%view(2)
axis tight
%caxis([0 1])
ylabel('Sorted Samples')
xlabel('Period [min]')

figReporter = Figure(fig);
figReporter.SnapshotFormat = 'tif';
figImg = Image(getSnapshotImage(figReporter,rpt));
figImg.Style = [figImg.Style {ScaleToFit}];
add(rpt,figImg)

close gcf 
    
close(rpt);

end

