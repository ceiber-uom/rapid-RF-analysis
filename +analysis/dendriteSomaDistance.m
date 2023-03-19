
function result = dendriteSomaDistance( filename, varargin )
% summary = extractDistanceMeasurements( filename, ... )
% 
% Derive relationship between "as the crow flies" distance and actual
%   (path-integral along the dendrite) distance (in 2D and 3D). 
% 
% The returned data structure has the following fields: 
%  .filename    - input filename
%  .soma        - selected soma location
%  .dendrite    - copy of the input anatomy XYZ dendrite locations
%  .distance_1d - geodesic (as-the-crow-flies) distance from the soma to
%                 each point in the dendritic tree (computed in 2D)
%  .distance_2d - path-length integral distance (in 2D) from the soma to
%                 each point in the dendritic tree
%  .distance_3d - path-length integral distance (in 3D) from the soma to
%                 each point in the dendritic tree. 
%  .stats       - Scholl analysis of dendritic tree 
%                 https://en.wikipedia.org/wiki/Sholl_analysis
% 
% Does not correctly handle cells with multle primary dendrites as built,
%  but I've implemented a workaround using -repeat [n] mode, which allows you
%  to pick a cell and analyse up to N primary dendrites. When you wish to
%  stop, push 'esc' when asked to select the next primary dendrite. 
%  Use 'tab' to swap the Y and Z axes while selecting a primary dendrite. 
% 
% This code is comfortable handling three kinds of input data: 
% - .mat files corresponding to the skeletonised retinal ganglion cells 
%        published in Bae et al. (2018)
%       (https://doi.org/10.1016/j.cell.2018.04.040)
%       (https://github.com/seung-lab/e2198-gc-analysis) 
% - .mat files traced in matlab, published by Eiber et al. (2019)
% - .hoc files exported from traced RGCs using Imaris 
%
% you can also pre-load an anatomy file using tools.loadAnatomy and call
%   result = analysis.dendriteSomaDistance( anat, ... )
% 
% Options
% -repeat [n] - automatically repeat the analysis N times (1 for each
%               primary dendrite) and merge the data 
% -voxel-size [92 66 66] nm - from Seung lab documentation, see loadAnatomy
%                             I think this is only applied to .mat files. 
% -scale [] - same as -voxel-size but in units of um. 
% -load-opts {} - additional options to tools.loadAnatomy()
% 
% -filter [3.5] - for summary fit line in red, ignore points further than
%                 this ratio (probably bad guess as to primary dendrite)
% -plot - make final output plot (enabled by default). You can also re-plot
%         the analysis results using analysis.DSD( result, '-plot' )
% -no-plot - suppress final figure generation
% -eab-fig - generate an additional figure (EAB September 2022)
% -debug - make progressive debug animation of distance calculation
% -soma [x y z] - set soma position initial guess
% -auto - automatically guess the soma position from minimum Z
% 
% v0.1 - 5 September 2022, Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

if nargin > 0 && strcmpi(filename,'all')
    loop_over_all_files(varargin{:})
    return
end

if nargin > 0 && any(named('-plot')) && isstruct(filename)
                            make_figure(filename), return
elseif any(named('-plot')), make_figure(get_('-plot')), return
end

figure(2)
if nargin == 0 || ischar(filename) && ~exist(filename,'file') 
    filename = pick_file;
end

if any(named('-repeat')) % repeat mode
    n_reps = get_('-repeat');
    varargin(find(named('-repeat')) + [0 1]) = []; % remove arg
    result = run_repeated_analysis(n_reps, filename, varargin{:});
    return
end

if any(named('-vo')), load_opts = {'-scale',get_('-vo')/1e3};
elseif any(named('-sc')), load_opts = {'-scale',get_('-sc')};
else load_opts = {}; 
end
if any(named('-load-o')), load_opts = [load_opts get_('-load-o')]; end

if isstruct(filename), s = filename; filename = s.name; 
else
    fprintf('Loading %s\n', filename )
    s = tools.loadAnatomy(filename, load_opts{:}); % s for skeleton
end

%%

seg_length_3D = sqrt( sum((s.node(s.edge(:,1),:) - ...
                           s.node(s.edge(:,2),:)).^2, 2)); 

seg_length_2D = sqrt( sum((s.node(s.edge(:,1),1:2) - ...
                           s.node(s.edge(:,2),1:2)).^2, 2)); 


if any(named('-soma')), soma_xy = get_('-soma'); 
    if numel(soma_xy) == 2, soma_xy(3) = min(s.node(:,3)); end
    [~,soma_id] = min(sum( (s.node - soma_xy).^2, 2));    
elseif isfield(s,'soma') && any(named('-auto'))
    [~,soma_id] = min(sum( (s.node - s.soma).^2, 2));
elseif any(named('-auto'))
    [~,soma_id] = min(s.node(:,3)); 
else
    %% Pick soma manually
    clf
    i_color = s.node(:,2);
    if any(named('-init')), i_color = get_('-init');
      if isempty(i_color),  i_color = s.node(:,2);      end
      if isstruct(i_color), i_color = i_color.distance_3d; end
    end

    seq = [1 3 2];
    if length(unique(s.node(:,3))) == 1, seq = [1 2 3]; end

    im = scatter(s.node(:,seq(1)), s.node(:,seq(2)), [], i_color, '.', ...
                     'UserData',s.node(:,seq(3)));
    axis image, grid on
    title('click to select, tab to swap Y/Z, esc to cancel')
    
    while true
        [x,y,b] = ginput(1);
        if isempty(b), result = []; return, end 
        % enter key allows select with more careful cursor
        if b == 9 % on tab key swap Y/Z
          ud = im.UserData;
          im.UserData = im.YData;
          im.YData = ud;
          im.CData = im.UserData;
          seq = seq([1 3 2]);
          continue  
        end
        if b > 3, result = []; return, end  % cancelled e.g. esc key
        break % any other key (i.e. mouse left/middle/right) is finish
    end

    soma_xy = [x y];
    [~,soma_id] = min(sum( (s.node(:,seq(1:2)) - soma_xy).^2, 2));        
    % soma_xy = [x y min(s.node(:,3))];
end


dist_to_soma_3D = inf * s.node(:,1); 
dist_to_soma_3D(soma_id) = 0;
dist_to_soma_2D = dist_to_soma_3D;
dist_to_soma_1D = sqrt(sum((s.node(:,1:2) - s.node(soma_id,1:2)).^2,2)); 

next = soma_id; 

do_debug_plot = any(named('-debug'));

if do_debug_plot
    %%
    clf
    h = scatter3(s.node(:,1), s.node(:,2), s.node(:,3), [], dist_to_soma_3D, '.');
    axis image, grid on
end

twocol_ = @(x) reshape(x,[],2); 

printInfo();

EAB_node_xyz = []; 

while any(~isfinite(dist_to_soma_3D))
    
    printInfo('Computing distance to soma [%0.1f%%] ... ', ...
                        100*mean(isfinite(dist_to_soma_3D)));
    
    if ~any(next)
        %% Need to jump over a gap
    
        sel = ~isfinite(dist_to_soma_3D); % what is not selected?
        
        inc_ids = find(~sel);
        exc_ids = find(sel); 
        
        [gap_dist,nnid] = arrayfun(@(u) min(sum((s.node(~sel,:) - ...
                                                 s.node(u,:)).^2, 2)), ...
                                                 exc_ids);
        
        [~,bnid] = min(gap_dist);
        
        nnid = inc_ids(nnid(bnid));
        gap_dist = sqrt(gap_dist(bnid));
        bnid = exc_ids(bnid); 
        
       
        gap_dist_2D = sqrt(sum((s.node(bnid,1:2)-s.node(nnid,1:2)).^2,2)); 
        
        dist_to_soma_3D(bnid) = dist_to_soma_3D(nnid) + gap_dist;
        dist_to_soma_2D(bnid) = dist_to_soma_2D(nnid) + gap_dist_2D;
        
        assert(isfinite(dist_to_soma_3D(bnid)))        
        next = bnid; 
        
        % sel = any(ismember(s.edge, bnid), 2);
        % node_ids = s.edge(sel,:);
        % next = node_ids( ~isfinite(dist_to_soma_3D(node_ids)));

    else        
        %% connection well defined for us already

        sel = any(ismember(s.edge, next), 2); % the edges which contain a node marked 'next'    
        node_ids = s.edge(sel,:); % the IDs of each node in one of the above edges
        
        if any(named('-eab-fig'))
          for rr = 1:size(node_ids,1)
            EAB_node_xyz = [EAB_node_xyz;
                            s.node(node_ids(rr,:),:); 
                            NaN NaN NaN]; %#ok<AGROW> 
          end
        end

        % next iteration, the list of nodes to analyse are the nodes in the
        % above collection of edges which weren't connected previously.
        % 
        % need to assign before updating distances because we'll forget which
        % ones weren't connected before
        next = node_ids( ~isfinite(dist_to_soma_3D(node_ids))); 

        % update the distances to the soma. the nodes in "node_ids" which have
        % already been assigned will be unchanged but the unconnected nodes
        % will be set to the distance of their connected neighbor  plus the
        % length of the connecting segment between the nodes    
        
        d = min( twocol_(dist_to_soma_3D(node_ids)), ...        
                                min( twocol_(dist_to_soma_3D(node_ids)), [], 2 ) + ... 
                                  seg_length_3D(sel) * [1 1]);
        assert(all(isfinite(d(:))))
        
        dist_to_soma_3D(node_ids) = d;

        dist_to_soma_2D(node_ids) = min( twocol_(dist_to_soma_2D(node_ids)), ...        
                                    min( twocol_(dist_to_soma_2D(node_ids)), [], 2 ) + ... 
                                  seg_length_2D(sel) * [1 1]);  
                              
    end
    
    if do_debug_plot
        h.CData = dist_to_soma_3D;
        % pause(0.01)
        caxis([0 max(dist_to_soma_3D(isfinite(dist_to_soma_3D)))])
    end
        
end

fprintf('Done! \n')
%%

% Putting in the following is a quick way to make a figure but it's not
% very good programming practice for building tools which are supposed to
% be general-purpose (i.e. use-many-times code as opposed to use-once code)
% 
% If you absolutely need to make figures like this, the thing you can do is
% put the figure in an "if any(named('-make-my-figure'))" block so that you
% can request the code to generate that figure but it doesn't break the
% code for other uses of the code, like below 
% 
% f2 = figure;
% plot3(pnodes(:,1),pnodes(:,2),pnodes(:,3),'LineWidth',2);
% axis image, grid on
% ca = gca;
% ca.XLim = [50,260];
% ca.YLim = [85,230];
% ca.XTick = 50:10:260;
% ca.YTick = 80:10:230;
% ln = ca.Children;
% xd = ln.XData;
% yd = ln.YData;
% 
% BECOMES (without the explicit XLim/YLim): 

if any(named('-eab-fig')) 
    figure;
    plot3(EAB_node_xyz(:,1),EAB_node_xyz(:,2),EAB_node_xyz(:,3),'LineWidth',2);
    axis image, grid on
    axis(axis * [1.05 -0.05 0 0; -0.05 1.05 0 0; ... % expand axes 5% 
                 0 0 1.05 -0.05; 0 0 -0.05 1.05]) 
    set(gca,'XTick', round(min(xlim),-1):10:max(xlim), ... % 10 um grid
            'YTick', round(min(ylim),-1):10:max(ylim))
end


result = struct; 
result.filename = filename; 
result.soma = s.node(soma_id,:); 
result.dendrite    = s.node; 
result.distance_1d = dist_to_soma_1D;
result.distance_2d = dist_to_soma_2D;
result.distance_3d = dist_to_soma_3D;
result.input_options = varargin;

cutoff = 3.5;
if any(named('-filt')), cutoff = get_('-filt'); end
result = compute_summary(result, cutoff, s); 

if ~any(named('-no-p')), make_figure(result, named), end
return




%% Scripts for looping analysis (-all and -repeat n)
function loop_over_all_files(varargin)

list = dir('../skeletons_GC/skel_*.mat'); 

warning(['Because of the primary dendrite problem, this code won''t ' ...
         'just work and give you all valid results'])

p_ = @(x) [x.folder filesep x.name];

summary = []; 

for ii = 1:numel(list)
    
    s = analysis.dendriteSomaDistance( p_(list(ii)), varargin{:}); 
    
    if isempty(summary), summary = s;
    else summary(end+1) = s;
    end
end

save('extractDistanceMeasurements_EB.mat','summary')

function result = run_repeated_analysis(n_replicates,varargin)

this = []; 
result = []; 

if ~exist(varargin{1},'file'), varargin{1} = pick_file; end

for iter = 1:n_replicates
    this = analysis.dendriteSomaDistance( varargin{:}, '-init', this );
    if isempty(this), break, end % ended by user

    pause(0.05)
    
    if isempty(result), result = this;
    else 
      for dist = {'distance_2d','distance_3d'}
        result.(dist{1}) = min(this.(dist{1}), result.(dist{1}));
      end
      result.distance_1d = result.distance_1d.*(iter-1)./(iter) + ...
                           this.distance_1d./(iter);       
    end 
end

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

cutoff = 3.5;
if any(named('-filt')), cutoff = get_('-filt'); end
result = compute_summary(result, cutoff); 

analysis.dendriteSomaDistance(result, '-plot')


%% Pick a file (peristent file path)
function filename = pick_file

  persistent f_path
  if isempty(f_path) || all(f_path == 0)
    f_path = '../skeletons_GC/';
  end

  [fn,fp] = uigetfile({'*.hoc';'*.mat'},'',f_path);
  if all(fn == 0), error('file selection cancelled'), end
  f_path = fp;
  filename = [fp fn]; 
return

%% Generate output figure (1D vs 2D, 3D metrics)
function make_figure(s, named)

if nargin < 2, named = @(x) false; end

xyz = s.dendrite;

C = lines(7);
clf

subplot(2,1,1)
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), [], s.distance_3d, '.');
axis image, tidyPlotForIllustrator, grid on %#ok<*DUALC>
ylabel(colorbar,'3D path distance (µm)')
% try , end, grid on

d_style = {'.','Color',[.3 .3 .3]};

subplot(2,3,4)    
plot(s.distance_1d, s.distance_3d, d_style{:})
xlabel('µm cartesan distance'), ylabel('µm path distance (3D)')
axis image, tidyPlotForIllustrator, grid on
set(gca,'userdata','1v3')
% try tidyPlotForIllustrator, end

if isfield(s,'stats')
    hold on
    errorbar( s.stats.x, s.stats.fit_1to3d_avg, ...
                         s.stats.fit_1to3d_std,'LineWidth',1.2)
    plot( s.stats.x, s.stats.x * s.stats.fit_1to3d_cutoff, ... 
                         '-','Color',[0 0 0 0.3])
end

if any(named('-2v3')) || ~isfield(s,'stats') % missing stats

    subplot(2,3,5)    
    plot(s.distance_1d, s.distance_2d, d_style{:})
    xlabel('µm cartesan distance'), ylabel('µm path distance (2D)')
    axis image, tidyPlotForIllustrator, grid on
    set(gca,'userdata','1v2')
    % try tidyPlotForIllustrator, end
    
    subplot(2,3,6)
    plot(s.distance_2d, s.distance_3d-s.distance_2d, d_style{:})
    xlabel('µm path distance (2D)'), ylabel('µm difference in distance (3D-2D)')
    axis image, tidyPlotForIllustrator, grid on
    set(gca,'userdata','2v3')
    return
else
    p = get(gca,'Position')./[1 1 1 4];
    set(gca,'Position', p.*[1 1 1 3]);
    axes('Position',p + [0 3*p(4) 0 0]);

    plot(s.distance_1d, s.distance_3d-s.distance_2d, d_style{:})
    ylabel('(3D-2D)')
    axis image, tidyPlotForIllustrator, grid on
    set(gca,'userdata','2v3')
end

%% Plot scholl distance
if isfield(s.stats, 'schollCount')
    %%
    subplot(2,3,5), cla

    smooth = @(y) conv(y([1 1:end end]), [1 1 1]/3,'valid');
    area(s.stats.schollRadius, smooth(s.stats.schollCount), ...
           'LineWidth',1.2,'EdgeColor',C(1,:),'FaceAlpha',0.3)
    xlabel('µm radius'), ylabel('Scholl Crossing Density')
    tidyPlotForIllustrator
    set(gca,'userdata','schollCount')

    p = get(gca,'Position')./[1 1 1 4];
    set(gca,'Position', p.*[1 1 1 3]);
    
    subplot(2,3,6), cla, hold on
    area(s.stats.schollRadius, smooth(s.stats.schollLength_1D), ...
           'LineWidth',1.2,'EdgeColor',C(1,:),'FaceAlpha',0.3)
    area(s.stats.schollRadius, smooth(s.stats.schollLength_3D), ...
           'LineWidth',1.2,'EdgeColor',C(2,:),'FaceAlpha',0.3)
    
    xlabel('µm radius'), ylabel('µm Scholl Length')
    tidyPlotForIllustrator
    set(gca,'userdata','schollCount')

    p = get(gca,'Position')./[1 1 1 4];
    set(gca,'Position', p.*[1 1 1 3]);
    legend('1D','3D','location','best'), legend boxoff
    
end


return

%% Compute summary metrics from data
function r = compute_summary(r, fit_cutoff, s)
% Compute 3D-to-1D relationship, density metrics, and Scholl analysis

fit_x = 0 : 5 : max(r.distance_1d);
fit_dx = mean(diff(fit_x)); 
fit_ok = (r.distance_3d < fit_cutoff * r.distance_1d); 
 
moving_fun = @(f,y) arrayfun(@(u) f(y( fit_ok & ...
                                abs(r.distance_1d-u) < fit_dx)), fit_x);

stats.x = fit_x;
stats.fit_1to3d_avg = moving_fun( @mean, r.distance_3d ); 
stats.fit_1to3d_std = moving_fun( @std, r.distance_3d ); 
stats.fit_1to3d_cutoff = fit_cutoff;
r.stats = stats;

if nargin < 3, return, end

%% 
% https://en.wikipedia.org/wiki/Sholl_analysis


radius = sqrt(sum((r.dendrite-r.soma).^2,2));
radius = radius(s.edge);

sr = max(max(r.distance_3d), max(radius(:)));
sr = 0: 1 : sr;

schollCount = @(r) sum( any(radius <= r, 2) & any(radius >  r, 2));
stats.schollRadius = sr;
stats.schollCount = arrayfun(schollCount,sr);

seg_length_3D = sqrt( sum((s.node(s.edge(:,1),:) - s.node(s.edge(:,2),:)).^2, 2)); 
node_len = 0*r.distance_1d;

for ii = 1:numel(seg_length_3D)
    node_len(s.edge(ii,:)) = node_len(s.edge(ii,:)) + seg_length_3D(ii)/2;
end

sum_in_ = @(n,x) sum(node_len(n <= x));
stats.schollLength_1D = arrayfun(@(x) sum_in_(r.distance_1d, x), sr );
stats.schollLength_3D = arrayfun(@(x) sum_in_(r.distance_3d, x), sr );

stats.schollLength_1D = diff([0 stats.schollLength_1D]);
stats.schollLength_3D = diff([0 stats.schollLength_3D]);

r.stats = stats; 

return

