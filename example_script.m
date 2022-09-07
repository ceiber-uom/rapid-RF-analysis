
% some of the best files:

% 20180703_Cell_1 #9
% 20180803_Cell_1 #1

clear
% d = utils.load('20180703_Cell_1 #9','-dir','../HEKA Radon/MAT','-pca');

d = utils.load('20210601_Cell_1 #12','-pca');

plots.standardFigure('Name','Standard PCA analysis'), clf
rdat = plots.plot_radon_IMG(d); 

% d = utils.prepareRadon(d, '-append'); 
% r = analysis.inverseRadon(d); 

%%
analysis.fitGaussianModel(d)




% plots.standardFigure('Name','Latency estimate'), clf
% t = analysis.estimateWaveLag(d.response_waves(:,1), d.time, d.expoData,'-plot'); 

%% Select time points for display in "total RF" figure

nK = size(d.response_waves,2); 
% selection based on maximum wave amplitude (in complex sense) of each
% PCA (or NNMF, if so inclined) component

% only consider maximal points in this window (could have edge effects)
window = (d.time > -0.1 & d.time < 0.9*max(d.time))'; 
[~,tt_points] = arrayfun(@(k) max(abs(hilbert(d.response_waves(:,k))) ... 
                                      .* window), 1:nK ); 
tt_points = sort([t.zero_index t.index tt_points]); 

%% Generate plot of 'total RF' at each time-point

plots.totalSensitivity(d, tt_points, '-im', rdat)
