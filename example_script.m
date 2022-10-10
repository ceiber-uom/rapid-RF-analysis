
% some of the best files:

% 20180703_Cell_1 #9
% 20180803_Cell_1 #1
% 20190904_Cell_02#16 OFF cell

clear
% d = tools.load('?','-nnmf', '-psth');
% d = tools.load('?','-pca');
d = tools.load('..\MAT\20220811_Cell_02#9[Radon_Flicker_ACH].mat','-PSTH','-pca');


plots.standardFigure('Name','Standard PCA analysis'), clf
rdat = plots.plot_radon_IMG(d); 

%%

% if ~exist('anat','var')
%     anat = tools.loadAnatomy();
% end
% analysis.dendriticDensity(anat, rdat,'-align');
% git
% d = tools.prepareRadon(d, '-append'); 
% r = analysis.inverseRadon(d); 

%% Demonstrate fitGaussianModel

plots.standardFigure('Name','Gaussian Model'), clf
gm = analysis.fitGaussianModel(d, '-nG',2,'-images','-use-c',[2 3]);

%% Demonstrate estimateWaveLag (not super reliable)

% plots.standardFigure('Name','Latency estimate'), clf
% t = analysis.estimateWaveLag(d.response_waves(:,1), d.time, d.expoData,'-plot'); 
% t = analysis.estimateWaveLag(d)

%% Demonstrate totalSensitivity

% plots.totalSensitivity(d, '-row',5,'-t', -0.1:0.1:0.8 )
plots.totalSensitivity(d, '-interactive' )

%% Demonstrate analysis.prediction (spots and annuli)

analysis.prediction(d, rdat)
