function [] = plotDifferencePSD(PSD,f,Pmax,n0,datasetName,sampPeriod,binSize,zeroPadding,derivativeIndicator)
%will plot difference PSD(1/f) spectra for pairs of datasets
% PSD: PSD amplitudes for each sample in each dataset
% f: Fourier frequencies for each sample in each dataset
% Pmax: period [min] of FFT peak with maximum PSD for each sample in each dataset
% n0: the number of dataset to be subtructed

outputPath = '/Users/muxika/MATLAB/result';
mkdir(outputPath);
outputPathSpectr = [outputPath '/differercePSD.pdf'];
Ndatasets = size(PSD,2);
S = cellfun(@(x) size(x,2), PSD);

%create report object
import mlreportgen.dom.*
import mlreportgen.report.*
rpt = Report(outputPathSpectr,'pdf'); 

%create table with FFT parameters for the whole dataset
parameters = cell(Ndatasets+1,5);
parameters(1,1:5) = {'dataset' 'binSize[mHz]' 'zero_padding' 'FFT_on_derivative' 'sampling_period[sec]'}; 
for p = 1:Ndatasets
    parameters(p+1,1) = {datasetName{p}};
    if size(binSize,2) ~= 1
        minBin = min(binSize);
        maxBin = max(binSize);
        binRange = {['bin size is in the range between ' num2str(minBin) ' and ' num2str(maxBin)]};
        parameters(p+1,2) = binRange;
    else
        parameters(p+1,2) = {num2str(binSize)};
    end
    parameters(p+1,3) = {num2str(zeroPadding)};
    parameters(p+1,4) = {num2str(derivativeIndicator)};
    parameters(p+1,5) = {num2str(sampPeriod(p))};
end
parameters(Ndatasets+2,1) = {'^{*} Samples in datasets are sorted by the period of the most intense PSD peak'};

table = Table(parameters);
add(rpt,table);
br = PageBreak();
add(rpt,br);



for n = 1:Ndatasets
        if n ~= n0
            fig = figure('visible','off');
            fig.Position = [1000 570 560 768];
            hold on;
            %cut bigger dataset
            nmin    = min([S(n0) S(n)]);
            f0      = f{n0}(1,1:nmin);            
            PSD0    = PSD{n0}(1,1:nmin);
            Pmax0   = Pmax{n0}(1,1:nmin);
            PSD1    = PSD{n}(1,1:nmin);
            Pmax1   = Pmax{n}(1,1:nmin);
            
            %sort the samples in dataset by the period of the most intense PSD peak
            [~,ind0]    = sort(Pmax0);
            [~,ind1]    = sort(Pmax1);            
            PSD0        = PSD0(ind0);
            PSD1        = PSD1(ind1);
            f0          = f0(ind0);
            %create a grid for surface plot
            [X,Y] = meshgrid((1000/60)./f0{1},1:nmin);
            %calculate difference PSD
            PSDdiff = cellfun(@minus,PSD1,PSD0,'UniformOutput',false);
            
            PSDmin = min(cellfun(@min, [PSD0 PSD1 PSDdiff]));
            PSDmax = max(cellfun(@max, [PSD0 PSD1 PSDdiff]));
            % plot PSD surface for dataset1
            ax1 = subplot(3,1,1);
            ax1.Position = [0.13,0.709,0.750,0.215];
            s = surf(X,Y,cell2mat(PSD1)');
            s.EdgeColor = 'none';
            c = colorbar;
            c.Label.String = 'PSD';
            %view(-18.698935578363781,72.538181818181812)
            view(2)
            axis tight
            if derivativeIndicator == 1
                xlim([1 14])
            end
            zlim([0 PSDmax])
            caxis([0 PSDmax])
            ylabel('Sorted Samples')
            xlabel('Period [min]')
            title(['1: ' strrep(datasetName{n},'_',' ')])
            
            % plot PSD surface for dataset0
            ax2 = subplot(3,1,2);
            ax2.Position = [0.13,0.409,0.750,0.215];
            s = surf(X,Y,cell2mat(PSD0)');
            s.EdgeColor = 'none';
            c = colorbar;
            c.Label.String = 'PSD';
            %view(-18.698935578363781,72.538181818181812)
            view(2)
            axis tight
            if derivativeIndicator == 1
                xlim([1 14])
            end
            zlim([0 PSDmax])
            caxis([0 PSDmax])
            ylabel('Sorted Samples')
            xlabel('Period [min]')
            title(['2: ' strrep(datasetName{n0},'_',' ')])
            
            % plot difference PSD surface for dataset1-dataset0
            ax3 = subplot(3,1,3);
            ax3.Position = [0.13,0.11,0.750,0.215];
            s = surf(X,Y,cell2mat(PSDdiff)');
            s.EdgeColor = 'none';
            c = colorbar;
            c.Label.String = 'Difference PSD';
            %view(-18.698935578363781,72.538181818181812)
            view(2)
            axis tight
            if derivativeIndicator == 1
                xlim([1 14])
            end
            zlim([PSDmin PSDmax])
            caxis([PSDmin PSDmax])
            ylabel('Sorted Samples')
            xlabel('Period [min]')
            title('1-2')
            
            figReporter = Figure(fig);
            figReporter.SnapshotFormat = 'tif';
            figImg = Image(getSnapshotImage(figReporter,rpt));
            figImg.Style = [figImg.Style {ScaleToFit}];
            add(rpt,figImg)

            close gcf 
        end   
end
    
close(rpt);

end

