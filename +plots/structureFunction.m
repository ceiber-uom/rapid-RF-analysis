
function structureFunction(anat, rdat, varargin)
% plots.structureFunction(anat, rdat, [anat_distance], ... )
% 
% Generate structure-function plot of dendrite density and distance
%   to soma vs receptive field strength. 
% 
% This code integrates the outputs of plots.anatomy,
%   analysis.dendriticDensity, and analysis.dendriteSomaDistance to
%   generate an overall results figure. The output figure has 4 panels:
% 
%   A - cell anatomy shown in context of receptive field (plots.anatomy)
%   B - anatomy and dendritic field density (analysis.dendriticDensity)
%   C - receptive field mapped onto the dendritic field points
%        (same colourscale as A)
%   D - scatterplot of 3D dendritic length (from analysis.dendriteSomaDist)
%        against receptive field strength, colour-coded by local dendritic
%        density
% 
%  The axis values for (D) can be customised as follows: 
%  -3d    : X-axis is 3D integrated dendrite path-length (default), 
%  -2d    : X-axis is 2D integrated dendrite path-length 
%  -1d    : X-axis is geodesic distance from dendrite point to soma
%  -depth : X-axis is the depth of the dendrite in the inner plexi layer
%  -local : X-axis is the local dendritic density. If this option is
%            selected, the color axis becomes whatever X would have been 
%           (i.e. the above options instead control the color axis). 
% 
% Other Options: 
% -dist [dm]    : use precomputed distance metrics
%                 (equivalent to third positional argument)
% -interactive  : Enable click-to-move cell
% -id [1]       : Select component ID to plot
% -z [0 0]      : set cell center (overrides -interactive)
% -xy [0 0]     : set cell center (works with -interactive)
% -anat-opts {} : pass additional options to plots.anatomy()
% -dd-opts   {} : pass additional options to analysis.dendriticDensity()
% -dsd-opts {}  : pass additional options to analysis.dendriteSomaDistance
% -radius [0]   : instead of taking the (interpolated) value of the RF at 
%                  each point in the dendritic tree, take the average 
%                  within some radius around each point.
% -d-blur [0]   : apply gaussian radius blur to computed dendritic
%                  density function (default: no blur).
% 
% v0.1 - 24 December 2022 - Calvin Eiber <c.eiber@ieee.org> 
% Changelog: Merry Christmas! 

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

xy_zero = [0 0];
if any(named('-z')), xy_zero = get_('-z'); 
elseif any(named('-xy')), xy_zero = get_('-xy'); 
end

anat_c = median(anat.node);

PA_opts = {}; 
if any(named('-anat-opts')), PA_opts = get_('-anat-opts'); end
if any(named('-id')), PA_opts = [PA_opts {get_('-id')}]; end

DD_opts = {}; 
if any(named('-dd-opts')), DD_opts = {get_('-dd-opts')}; end
%%

clf
cf = gcf;
cf.Position = [-1026,390,981,1112];

subplot(3,2,1)
% plots.anatomy(anat, rdat, PA_opts{:}, '-xy', xy_zero, '-clf')
plots.anatomy(anat, rdat, PA_opts{:}, '-clf')


h = findobj(gca,'type','image');
[gx,gy] = ndgrid(h(1).XData, h(1).YData);
% afi = griddedInterpolant(gx, gy, h(1).CData);

subplot(3,2,2)
analysis.dendriticDensity(anat, rdat, DD_opts{:}, '-z', xy_zero-anat_c(1:2))
d = findobj(gca,'type','image'); 

%%

if any(named('-dist')), metrics = get_('-dist'); 
elseif nargin>2 && isstruct(varargin{1}), metrics = varargin{1};
else
    %%
    DSD_opts = {}; 
    if any(named('-dsd-opts')), DSD_opts = get_('-dsd-opts'); end
    if any(named('-repeat')), DSD_opts = [DSD_opts {get_('-repeat')}]; end

    f = figure; 
    metrics = analysis.dendriteSomaDistance(anat.name, DSD_opts{:}, '-no-plot');
    close(f)
end

%% Get RF value near each node

radius = []; 
if any(named('-r')), radius = get_('-r'); end
anat.xyz = anat.node - anat_c + [xy_zero 0]; 

if isempty(radius)
  
  % get directly the value of the RF at each dendrite point
  afi = griddedInterpolant(gx, gy, h(1).CData);
  dendrite_field_value = afi( anat.xyz(:,1), anat.xyz(:,2));

else 
  
  % At each dendrite point take the average in a radius around the dendrite
  radius_sq = radius.^2; 
  dendrite_field_value = zeros(size(anat.xyz(:,1)));
  rf_img = h(1).CData; 

  for pp = 1:size(anat.node,1)
    % distance to point squared
    d2p_sq = (gx-anat.xyz(pp,1)).^2 + (gy-anat.xyz(pp,2)).^2;
    sel = d2p_sq <= radius_sq; 
    if ~any(sel), [~,sel] = min(abs(d2p_sq(:))); end
    dendrite_field_value(pp) = mean(rf_img(sel)); 
  end
end

% get directly the value of the RF at each dendrite point

if any(named('-d-blur'))
    %% Optionally smooth the density image out 
    dbgr = get_('-d-blur');

    gauss = exp(- (gx.^2 + gy.^2)/dbgr.^2/2 );
    gauss = gauss/sum(gauss(:)); 

    im = conv2(d(1).CData, gauss,'same'); 
    d(1).CData = im; 
end

ddi = griddedInterpolant(gx, gy, d(1).CData');
metrics.local_density = ddi( anat.xyz(:,1), anat.xyz(:,2));

subplot(3,2,1)
disp('Select roi for analysis');
roi = drawrectangle;
pts = roi.Position;
rows = rdat.range > pts(1) & rdat.range < pts(3);
col = rdat.range > pts(2) & rdat.range < pts(4);
rdat_roi = rdat.images{get_('-id')}(rows,col);
dd_roi = d.CData(rows,col);

subplot(3,2,4)
scatter(dd_roi,rdat_roi,'b.')
xlabel('dendrite density');
ylabel('receptive field strength')



%% Show heatmap on the dendrites themselves (zoomed in)

subplot(3,2,5)
scatter(anat.node(:,1), anat.node(:,2), 3,  dendrite_field_value, 'filled')
axis image xy, hold on, caxis(h(1).Parent.CLim)
axis image xy off
ca = gca;
plot([ca.XLim(1),ca.XLim(1)+50],[ca.YLim(1)-10,ca.YLim(1)-10],'k-')
text(ca.XLim(1)+15,ca.YLim(1)-15,sprintf('50 %cm',char(181)));

try tidyPlotForIllustrator, end

%%

subplot(3,2,6)

x_value = 'distance_3d'; 

if any(named('-3d')), x_value = 'distance_3d';
elseif any(named('-2d')), x_value = 'distance_2d';
elseif any(named('-1d')), x_value = 'distance_1d'; 
elseif any(named('-depth')), x_value = 'depth';
    metrics.depth = metrics.dendrite(:,3); 
end

c_value = 'local_density';
if any(named('-loc')), c_value = x_value; 
    x_value = 'local_density';
end

scatter(metrics.(x_value), dendrite_field_value,3, metrics.(c_value),'filled')
% set(gca,'Position',get(gca,'Position') + [0.05 0 -0.1 0])
xl = [strrep(x_value,'_',' '),' (',char(181),'m)'];
xlabel(xl), ylabel('receptive field strength'); 
try tidyPlotForIllustrator, end

cb = colorbar; 
yl = sprintf('%s: SA (%cm^2/pixel) or V (%cm^3/pixel)',strrep(c_value,'_',' '),...
    char(181),char(181));
ylabel(cb,yl);

if ~any(named('-loc')), caxis(d(1).Parent.CLim); end

o = flipud(get(gcf,'Children')); 
cb.Position = [0.93 0.11 0.015 0.8];
o(2).Position([1 3]) = [0.07 0.015];
ylabel(o(2),'Receptive Field Strength')
set(gcf,'Name','Structure-Function Relationship')

if any(named('-interactive'))
    set(o(1).Children,'HitTest','off')
    o(1).ButtonDownFcn = @(a,b) on_axes_click(a,b, ...
                                    anat,rdat,metrics,varargin{:} );
    title(o(1),'click to move cell')
end

%%

return





function on_axes_click(~, event, varargin )
% regenerate the figure with a new XY point

xy = event.IntersectionPoint(1:2); 
named = @(n) strncmpi(varargin,n,length(n));

if any(named('-xy')), varargin(find(named('-xy'))+1) = {xy}; 
else varargin = [varargin {'-xy', xy}]; 
end

plots.structureFunction(varargin{:})

return
