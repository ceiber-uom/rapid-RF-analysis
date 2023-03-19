
% some of the best files:

% 20180703_Cell_1 #9
% 20180803_Cell_1 #1
% 20190904_Cell_02#16 OFF cell

clear
close all
% d = tools.load('?','-nnmf', '-psth');
% d = tools.load('?','-pca');
cell = '20211122_Cell_03#18[Radon_Flicker_ACH].mat';
p = ['..\MAT\',cell];


% d = tools.load(p,'-pca','-nK',3);
d = tools.load(p,'-PSTH','-pca','-nK',3);
% d = tools.load(p,'-PSTH','-pca','-nK',3,'-smooth');
% 20190904_Cell_02: '-smooth'

plots.standardFigure('Name','Standard PCA analysis'), clf
rdat = plots.plot_radon_IMG(d,'-units');
% Use flag '-units' for fitGaussianModel

%%

% Examples for dendriticDensity:
% 20190904_Cell_02#17
% close(1)
if ~exist('anat','var')
    anat = tools.loadAnatomy();
end
% 
% fi = 'U:\PRJ-vnrg2019\V19_PAPERS\V19_Elissa_Radon\IMARIS_RECON\20200116_Cell_02\20190904_Cell_02.hoc';
% result = analysis.dendriteSomaDistance(fi,'-repeat',3);
% result = analysis.dendriteSomaDistance(fi);

% di = analysis.dendriticDensity(anat, rdat,'-area-density');

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

% plots.standardFigure('Name','Gaussian Model'), clf
% gm = analysis.fitGaussianModel(d, '-nG',2,'-ortho');
% gm = analysis.fitGaussianModel(d,'-nG',2,'-images','-use-c',1:2,'-el');

% Convert amplitude into imp/s/pixel
% nP = length(gm.fit_params);
% nK = 3;
% idx = nP-nK+1:nP;
% for kk = 1:nK       
%     bs = mean(rdat.wave(rdat.time<0,kk),1);
%     dif = arrayfun(@(r) diff([rdat.wave(r,kk),bs]), 1:size(rdat.wave,1));
%     dif = dif';
%     [mx,imax] = max(abs(dif)); % mx: max increase or decrease from baseline   
%     gm.fit_params(:,idx(kk)) = gm.fit_params(:,idx(kk)).*mx;
% end

%% For presentation: Zoom in on RF map
% 20200116_Cell_02#15
% fn = 'Cell_Images\20200116_Cell_02.png';
% im = '20200116_Cell_02.png';
% date = str2double(im(1:8));
% T = readtable('U:\PRJ-VisionLab\Elissa_Belluccini\Spreadsheets\Imaris_Recon_Coordinates.xlsx');
% sel = strcmp(T.File,im);
% midx = T.cx(sel);
% midy = T.cy(sel);
% [anatomy,xy] = loadAnatomy_EB(fn,midx,midy);
% plots.plot_anatomy(date,rdat,anatomy,xy,gm,midx,midy)

%% Demonstrate estimateWaveLag (not super reliable)

% plots.standardFigure('Name','Latency estimate'), clf
% t = analysis.estimateWaveLag(d.response_waves(:,1), d.time, d.expoData,'-plot'); 
% t = analysis.estimateWaveLag(d)

%% Demonstrate plots.structureFunction
% 20190904_Cell_02#16 OFF'-repeat',3
% 20200116_Cell_02#15, ON '-repeat',2
% 20210601_Cell_01#12 ONOFF*
% 20210608_Cell_05#10 OFF '-Vm','-id',2 (components are different...)
% 20211122_Cell_02#9 ONOFF '-repeat',2
% 20211122_Cell_03#18 OFF
% 20211122_Cell_04#9 ONOFF '-repeat',3
% 20211129_Cell_02#14 OFF, '-repeat',5
% 20220228_Cell_02#9 OFF, '-repeat',3
% 20220401_Cell_02#8 OFF, '-repeat',2
% 20220524_Cell_01#14 OFF, '-Vm','-repeat',4
% 20220624_Cell_01#10 ON
% 20220624_Cell_03#6 OFF, '-repeat',2
% 20220825_Cell_03#6[Radon_Flicker_ACH_long] OFF, '-repeat',3
% 20220908_Cell_02#6[Radon_Flicker_ACH_long] OFF,'-nnmf', '-repeat',5
%
% plots.structureFunction(anat,rdat,'-anat-opts','-id',1,...
%     '-dd-opts','-area-density','-dsd-opts','-repeat',3,'-radius',24)
% plots.structureFunction(anat,rdat,'-anat-opts','-id',1,...
%     '-dd-opts','-area-density')
plots.structureFunction(anat,rdat,'-anat-opts','-id',1,...
    '-dd-opts','-area-density')%,'-dsd-opts','-repeat',2)

%% Demonstrate analysis.prediction (spots and annuli)

analysis.prediction(d, rdat)

%%

gm = analysis.fitGaussianModel(d,'-nG',2,'-images','-el','-use-c',1:3);
