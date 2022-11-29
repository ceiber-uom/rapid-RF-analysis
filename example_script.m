
% some of the best files:

% 20180703_Cell_1 #9
% 20180803_Cell_1 #1
% 20190904_Cell_02#16 OFF cell

clear
% d = tools.load('?','-nnmf', '-psth');
% d = tools.load('?','-pca');
cell = '20220531_Cell_02#2[Radon_Flicker_ACH].mat';
p = ['..\MAT\',cell];

d = tools.load(p,'-PSTH','-pca','-nK',3);

plots.standardFigure('Name','Standard PCA analysis'), clf
rdat = plots.plot_radon_IMG(d,'-units');
% Use flag '-units' for fitGaussianModel

%%

% if ~exist('anat','var')
%     anat = tools.loadAnatomy();
% end
% analysis.dendriticDensity(anat, rdat,'-align');
% git
% d = tools.prepareRadon(d, '-append'); 
% r = analysis.inverseRadon(d); 

%% Demonstrate totalSensitivity
% Examples for totalSensitivity:
% 20220811_Cell_02#9[Radon_Flicker_ACH] ON sustained
% 20220513_Cell_02#9
% 20211129_Cell_02#14[Radon_Flicker_ACH] OFF cell, black bar

% plots.totalSensitivity(d, '-row',5,'-t', -0.1:0.1:0.8 )
% plots.totalSensitivity(d, '-interactive' )

%% Demonstrate fitGaussianModel
% Examples for fitGaussianModel:
% 20220513_Cell_02#9[Radon_Flicker_ACH] ON sustained. 3 significant Gaussian fits
% 20220531_Cell_02#2
% 20220728_Cell_01#8 ONOFF cell. 

plots.standardFigure('Name','Gaussian Model'), clf
% gm = analysis.fitGaussianModel(d, '-nG',2,'-ortho');
gm = analysis.fitGaussianModel(d,'-nG',2,'-images','-use-c',1:3);

% Convert amplitude into imp/s/pixel
nP = length(gm.fit_params);
nK = 3;
idx = nP-nK+1:nP;
for kk = 1:nK       
    bs = mean(rdat.wave(rdat.time<0,kk),1);
    dif = arrayfun(@(r) diff([rdat.wave(r,kk),bs]), 1:size(rdat.wave,1));
    dif = dif';
    [mx,imax] = max(abs(dif)); % mx: max increase or decrease from baseline   
    gm.fit_params(:,idx(kk)) = gm.fit_params(:,idx(kk)).*mx;
end

%% For presentation: Zoom in on RF map
% 20200116_Cell_02#15
fn = 'Cell_Images\20200116_Cell_02.png';
im = '20200116_Cell_02.png';
date = str2double(im(1:8));
T = readtable('U:\PRJ-VisionLab\Elissa_Belluccini\Spreadsheets\Imaris_Recon_Coordinates.xlsx');
sel = strcmp(T.File,im);
midx = T.cx(sel);
midy = T.cy(sel);
[anatomy,xy] = loadAnatomy_EB(fn,midx,midy);
plots.plot_anatomy(date,rdat,anatomy,xy,gm,midx,midy)

%% Demonstrate estimateWaveLag (not super reliable)

% plots.standardFigure('Name','Latency estimate'), clf
% t = analysis.estimateWaveLag(d.response_waves(:,1), d.time, d.expoData,'-plot'); 
% t = analysis.estimateWaveLag(d)


%% Demonstrate analysis.prediction (spots and annuli)

analysis.prediction(d, rdat)
