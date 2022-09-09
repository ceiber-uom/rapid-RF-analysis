
% some of the best files:

% 20180703_Cell_1 #9
% 20180803_Cell_1 #1
% 20190904_Cell_02#16 OFF cell

clear
% d = utils.load('?','-nnmf', '-psth');
d = utils.load('?','-pca');

plots.standardFigure('Name','Standard PCA analysis'), clf
rdat = plots.plot_radon_IMG(d); 

% d = utils.prepareRadon(d, '-append'); 
% r = analysis.inverseRadon(d); 

%%
plots.standardFigure('Name','Gaussian Model'), clf
gm = analysis.fitGaussianModel(d, '-nG',2);

gm = gm(end); 

% plots.standardFigure('Name','Latency estimate'), clf
% t = analysis.estimateWaveLag(d.response_waves(:,1), d.time, d.expoData,'-plot'); 
% t = analysis.estimateWaveLag(d)

fitted_model = rdat;
fitted_model.y_all = gm.predicted_RF;

plots.standardFigure('Name','Gaussian Model output')
plots.plot_radon_IMG(fitted_model)

%% Add radii to radon image
h = get(gcf,'Children'); 
h = findobj(h,'YLim',h(1).YLim);

FWHM_circle = [cos(linspace(0,2*pi,61));
               sin(linspace(0,2*pi,61))]' * sqrt(2*log(2));
C = lines(7);                

for ii = 1:numel(h) % for each axis
    
    hold(h(ii),'on')
    for gg = 1:gm.n_gaussians

        xy = gm.center_xy(gg,:) + gm.guass_radius(gg) * FWHM_circle;
        plot(xy(:,2),xy(:,1),'-','Color',C(gg,:),'LineWidth',1.2,'Parent',h(ii))
        plot(gm.center_xy(gg,2),gm.center_xy(gg,1),'.','Color',C(gg,:),'Parent',h(ii))
    end
end
   

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

% plots.totalSensitivity(d, '-row',5,'-t', -0.1:0.1:0.8 )

plots.totalSensitivity(d, '-interactive' )