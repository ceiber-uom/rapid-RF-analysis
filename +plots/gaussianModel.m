

function gaussianModel(rdat, gm, varargin)
% plots.gaussianModel(radon data, gaussian model, ...)
% 
% View output 

if nargin == 1, gm = analysis.fitGaussianModel(rdat); end
if numel(gm) > 1, gm = gm(end); end

named = @(n) strncmpi(varargin,n,length(n));

rdat.y_all = gm.predicted_RF;

%% Add radii to radon image

clf
d = plots.plot_radon_IMG(rdat, varargin{:});
h = get(gcf,'Children'); 
h = findobj(h,'YLim',h(1).YLim); % select from children of gcf the radon images

theta = linspace(0,2*pi,61);
oneSD_circle = [cos(theta); sin(theta)]';
% FWHM_circle = that * sqrt(2*log(2));

C = lines(max([gm.n_gaussians]));
for ii = 1:numel(h) % for each axis
  
  hold(h(ii),'on')
  for gg = 1:gm.n_gaussians

    style = {'Color',C(gg,:),'LineWidth',1.2,'Parent',h(ii)};

    if isfield(gm,'gauss_eccentricity')
      ecc = sqrt(1 - (gm.gauss_eccentricity(gg) .* ...
                  cos(theta-deg2rad(gm.gauss_angle(gg))).^2));
      oneSD_circle = [cos(theta)./ecc; sin(theta)./ecc]';
    end

    xy = gm.center_xy(gg,:) + gm.gauss_radius(gg) * oneSD_circle;
    plot(xy(:,1),xy(:,2),'-',style{:})
    plot(gm.center_xy(gg,1),gm.center_xy(gg,2),'.',style{:})
  end
  axis(h(ii),'image','xy')
end

% Adjust also color limits to be consistent 
h = findobj(gcf,'type','axes');
h = h(cellfun(@(c) any(c~=[0 1]), {h.CLim}));
set(h,'CLim',[-1 1]*max(abs([h.CLim])));

% Move all the axes down to make room for timeseries trace
if any(named('-no-k')) || ~isfield(gm,'kinetic'), return, end

h = findobj(gcf,'type','axes');

for ii = 1:numel(h)
  h(ii).Position([2 4]) = h(ii).Position([2 4]) .* 0.75;
end

p = cat(1,h.Position); % axis positions 
p = [min(p(:,1)) max(p(:,[2 4])) * [1; 1.1] ...
     max(p(:,[1 3])*[1;1]) - min(p(:,1)) max(p(:,4))];
axes('Position', p)

plot(d.time, gm.kinetic,'LineWidth',1.2), hold on
if ~any(named('-no-b')), 
    plot(d.time, gm.resting, 'Color', [.5 .5 .5],'LineWidth',1.2);
end
plot(d.time([1 end]),[0 0],'Color',[0 0 0 0.3]);
tidyPlotForIllustrator, xlim(d.time([1 end]))
for ss = 1:rdat.nStimuli % add stim bars to the waves plot
    rectangle('Position',rdat.stim_bar(ss,0.1), ... 
              'FaceColor',[0 0 0 0.3], 'EdgeColor','none')
end


return

