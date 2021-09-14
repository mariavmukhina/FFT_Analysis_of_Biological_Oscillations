function [] = plot2DEnergyHistogram(E,P)
%plot 2D histogram of energies (powers) E of intense FFT peaks found by findpeaks
%function applied to PSD spectra of the samples in dataset taken under the
%same conditions; P is corresponding periods of oscillation

%find the cell with biggest number of elements
N = max(cellfun(@(c) size(c,1), E));
%pad remaining cells with corresponding number of NaNs
E = cell2mat(cellfun(@(c) [c;nan(N-size(c,1),1)],E,'UniformOutput', false));
P = cellfun(@(c) c',P,'UniformOutput', false);
P = cell2mat(cellfun(@(c) [c;nan(N-size(c,1),1)],P,'UniformOutput', false));

%prepare vectors for histogram
P = reshape(P, [size(P,1)*size(P,2),1]);
P(isnan(P)) = [];
E = reshape(E, [size(E,1)*size(E,2),1]);
E(isnan(E)) = [];

figure;
histogram2(P,E,'Normalization','probability','FaceColor','flat','DisplayStyle','tile');
colorbar
xlabel('Period of oscillations [min]')
ylabel('Energy')

end

