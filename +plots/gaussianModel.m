

function gaussianModel(rdat, gm, varargin)
% plots.gaussianModel(radon data, gaussian model, ...)
% 
% View output 

if nargin == 1, gm = analysis.fitGaussianModel(rdat); end
if numel(gm) > 1, gm = gm(end); end

rdat.y_all = gm.predicted_RF;

%% Add radii to radon image

clf
plots.plot_radon_IMG(rdat)
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
      ecc = sqrt(1 - (gm.gauss_eccentricity .* ...
                  cos(theta-deg2rad(gm.gauss_angle)).^2));
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

return

