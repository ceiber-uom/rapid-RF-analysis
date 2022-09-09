
function result = dendriteSomaDistance( filename, varargin )
% summary = extractDistanceMeasurements( filename, ... )
% 
% Derive relationship between "as the crow flies" distance and actual
%   (path-integral along the dendrite) distance (in 2D and 3D). 
% 
% Does not correctly handle cells with multle primary dendrites as built,
% but I've implemented a workaround using -rep [n] mode, which allows you
% to pick a cell and analyse up to N primary dendrites (when you wish to
% stop, push 'esc' or any keyboard key when selecting the next primary
% dendrite)
% 
% This code is comfortable handling three kinds of input data: 
% - .mat files corresponding to the skeletonised retinal ganglion cells 
%        published in Bae et al. (2018)
%       (https://doi.org/10.1016/j.cell.2018.04.040)
%       (https://github.com/seung-lab/e2198-gc-analysis) 
% - .mat files traced in matlab, published by Eiber et al. (2019)
% - .hoc files exported from traced RGCs using Imaris 
% 
% Options
% -repeat [n] - automatically repeat the analysis N times (1 for each
%               primary dendrite) and merge the data 
% -voxel-size [92 66 66] nm - from documentation [of what?]
% -filter [3.5] - for summary fit line in red, ignore points further than
%                 this ratio (probably bad guess as to primary dendrite)
% -plot - make final output plot (enabled if nargout == 0)
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

if nargin == 0 || ~exist(filename,'file'), filename = pick_file; end

if any(named('-re')) % repeat mode
    n_reps = get_('-re');
    varargin(find(named('-re')) + [0 1]) = []; % remove arg
    result = run_repeated_analysys(n_reps, filename, varargin{:});
    return
end

   

fprintf('Loading %s\n', filename )
s = get_data(filename); % s for skeleton

vox_to_um = [92 66 66] / 1e3; % from the readme
if isfield(s,'scale'), vox_to_um = s.scale; end
if any(named('-vo')),  vox_to_um = get_('-vo') / 1e3; end

s.n = s.n .* vox_to_um; % voxel coordinates of centerline tracing

%%

seg_length_3D = sqrt( sum((s.n(s.e(:,1),:) - s.n(s.e(:,2),:)).^2, 2)); 
seg_length_2D = sqrt( sum((s.n(s.e(:,1),1:2) - s.n(s.e(:,2),1:2)).^2, 2)); 


if any(named('-soma')), soma_xy = get_('-soma'); 
    if numel(soma_xy) == 2, soma_xy(3) = min(s.n(:,3)); end
    [~,soma_id] = min(sum( (s.n - soma_xy).^2, 2));    
elseif any(named('-auto')), [~,soma_id] = min(s.n(:,3)); 
else
    %% Pick soma manuyally
    clf
    i_color = s.n(:,2);
    if any(named('-init')), i_color = get_('-init');
      if isempty(i_color),  i_color = s.n(:,2);      end
      if isstruct(i_color), i_color = i_color.distance_3d; end
    end

    seq = [1 3 2];

    if length(unique(s.n(:,3))) == 1, seq = [1 2 3]; end

    im = scatter(s.n(:,seq(1)), s.n(:,seq(2)), [], i_color, '.', ...
                     'UserData',s.n(:,seq(3)));
    axis image, grid on
    
    b = 9;
    while b==9 % TAB key
        [x,y,b] = ginput(1);
        if isempty(b), result = []; return, end % enter key
        if b == 9 % on tab key swap Y/Z
          ud = im.UserData;
          im.UserData = im.YData;
          im.YData = ud;
          im.CData = im.UserData;
          seq = seq([1 3 2]);
          continue  
        end
        if b > 3, result = []; return, end  % cancelled e.g. esc key
        break
    end

    soma_xy = [x y];
    [~,soma_id] = min(sum( (s.n(:,seq(1:2)) - soma_xy).^2, 2));        
    % soma_xy = [x y min(s.n(:,3))];
end


dist_to_soma_3D = inf * s.n(:,1); 
dist_to_soma_3D(soma_id) = 0;
dist_to_soma_2D = dist_to_soma_3D;
dist_to_soma_1D = sqrt(sum((s.n(:,1:2) - s.n(soma_id,1:2)).^2,2)); 

next = soma_id; 

do_debug_plot = any(named('-debug'));

if do_debug_plot
    %%
    clf
    h = scatter3(s.n(:,1), s.n(:,2), s.n(:,3), [], dist_to_soma_3D, '.');
    axis image, grid on
end

twocol_ = @(x) reshape(x,[],2); 

printInfo();

while any(~isfinite(dist_to_soma_3D))
    
    printInfo('Computing distance to soma [%0.1f%%] ... ', ...
                        100*mean(isfinite(dist_to_soma_3D)));
    
    if ~any(next)
        %% Need to jump over a gap
    
        sel = ~isfinite(dist_to_soma_3D); % what is not selected?
        
        inc_ids = find(~sel);
        exc_ids = find(sel); 
        
        [gap_dist,nnid] = arrayfun(@(u) min(sum((s.n(~sel,:) - ...
                                                 s.n(u,:)).^2, 2)), ...
                                                 exc_ids);
        
        [~,bnid] = min(gap_dist);
        
        nnid = inc_ids(nnid(bnid));
        gap_dist = sqrt(gap_dist(bnid));
        bnid = exc_ids(bnid); 
        
       
        gap_dist_2D = sqrt(sum((s.n(bnid,1:2)-s.n(nnid,1:2)).^2,2)); 
        
        dist_to_soma_3D(bnid) = dist_to_soma_3D(nnid) + gap_dist;
        dist_to_soma_2D(bnid) = dist_to_soma_2D(nnid) + gap_dist_2D;
        
        assert(isfinite(dist_to_soma_3D(bnid)))        
        next = bnid; 
        
        % sel = any(ismember(s.e, bnid), 2);
        % node_ids = s.e(sel,:);
        % next = node_ids( ~isfinite(dist_to_soma_3D(node_ids)));

    else        
        %% connection well defined for us already

        sel = any(ismember(s.e, next), 2); % the edges which contain a node marked 'next'    
        node_ids = s.e(sel,:); % the IDs of each node in one of the above edges

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

result = struct; 

result.filename = filename; 

result.soma = s.n(soma_id,:); 
result.dendrite    = s.n; 
result.distance_1d = dist_to_soma_1D;
result.distance_2d = dist_to_soma_2D;
result.distance_3d = dist_to_soma_3D;
result.input_options = varargin;

if ~any(named('-no-p')), make_figure(result), end


%% Compute summary

cutoff = 3.5;
if any(named('-filt')), cutoff = get_('-filt'); end
result = compute_summary(result, cutoff); 

return


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

function result = run_repeated_analysys(n_replicates,varargin)

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

analysis.dendriteSomaDistance(result, '-plot')


named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

cutoff = 3.5;
if any(named('-filt')), cutoff = get_('-filt'); end
result = compute_summary(result, fit_cutoff); 


%% Pick a file (peristent file path)
function filename = pick_file

  persistent f_path
  if isempty(f_path) || all(f_path == 0)
    f_path = '../skeletons_GC/';
  end

  [fn,fp] = uigetfile({'skel*.mat';'anat*.mat';'*.hoc'},'',f_path);
  if all(fn == 0), error('file selection cancelled'), end
  f_path = fp;
  filename = [fp fn]; 
return

%% Generate output figure (1D vs 2D, 3D metrics)
function make_figure(s)

xyz = s.dendrite;

clf
subplot(2,1,1)
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), [], s.distance_3d, '.');
axis image, tidyPlotForIllustrator, grid on %#ok<*DUALC>
ylabel(colorbar,'3D path distance (�m)')
% try , end, grid on

subplot(2,3,4)    
plot(s.distance_1d, s.distance_3d,'.')
xlabel('�m cartesan distance'), ylabel('�m path distance (3D)')
axis image, tidyPlotForIllustrator, grid on
% try tidyPlotForIllustrator, end

subplot(2,3,5)    
plot(s.distance_1d, s.distance_2d,'.')
xlabel('�m cartesan distance'), ylabel('�m path distance (2D)')
axis image, tidyPlotForIllustrator, grid on
% try tidyPlotForIllustrator, end

subplot(2,3,6)
plot(s.distance_2d, s.distance_3d-s.distance_2d,'.')
xlabel('�m path distance (2D)'), ylabel('�m difference in distance (3D-2D)')
axis image, tidyPlotForIllustrator, grid on
% try tidyPlotForIllustrator, end

%% Compute summary metrics from data
function s = compute_summary(s, fit_cutoff)

fit_x = 0 : 5 : max(s.distance_1d);
fit_dx = mean(diff(fit_x)); 
fit_ok = (s.distance_3d < fit_cutoff * s.distance_1d); 
 
moving_fun = @(f,y) arrayfun(@(u) f(y( fit_ok & ...
                                abs(s.distance_1d-u) < fit_dx)), fit_x);

s.fit_1d = fit_x;
s.fit_3d_mean = moving_fun( @mean, s.distance_3d ); 
s.fit_3d_std = moving_fun( @std, s.distance_3d ); 

subplot(2,3,4), hold on
errorbar( fit_x, s.fit_3d_mean, s.fit_3d_std,'LineWidth',1.5)
plot( fit_x, fit_cutoff * fit_x, '-','Color',[0 0 0 0.3])




function skel = get_data(filename)

% TODO - determine how to correctly parse file
[~,stub,ext] = fileparts(filename);

if strcmp(ext,'.mat')

    vars = whos('-file',filename);
    
    if any(strcmp({vars.name},'dendrite'))

        disp(['Loading ' stub ext])
        error TODO_load_old_anatomy
        old = load(filename,'dendrite','soma');

        skel.n = []; 
        skel.e = []; 
        skel.f = []; 
        skel.scale = 1;

        % if isfield(old,'soma')

        

    elseif all(ismember({vars.name},'root','n','e'))
        disp(['Loading ' stub ext])
        skel = load(filename,'n','e'); 

    else error('unknown internal format: %s.mat', stub)
    end




elseif strcmp(ext,'.hoc') || strcmp(ext,'.geo')
  %% Parse HOC file
  disp(['Reading ' stub ext])
  fi = fopen(filename,'rt');
  oc = onCleanup(@() fclose(fi));

  skel = struct;
  skel.n = []; % node XYZ
  skel.e = []; % skeleton edges in graph
  skel.f = []; % filament ID
  skel.d = []; % filament diameter (not used)
  skel.scale = 1;

  obj_name = {};
  obj_idx = [];

  while ~feof(fi)
    s = fgetl(fi);
    s = regexprep(s,'//.*',''); % remove comments
    if isempty(s), continue, end
    
    if contains(s,'{')
        obj_name(end+1) = strtrim(regexp(s,'^[^{]*','match'));
        obj_idx = numel(obj_name);
        continue
    end
    
    if contains(s,'connect')
        val = str2double(regexp(s,'(?<=\()\d','match'));
        oid = find(cellfun(@(n) contains(s,n), obj_name));
        seq = strfind(s,obj_name{oid(1)}) > strfind(s,obj_name{oid(2)});
        if seq, oid = oid([2 1]); end

        if val(1), val(1) = find(skel.f == oid(1),1,'last');
        else       val(1) = find(skel.f == oid(1),1);
        end
        if val(2), val(2) = find(skel.f == oid(2),1,'last');
        else       val(2) = find(skel.f == oid(2),1);
        end
        skel.e = [skel.e; val];
        continue
    end

    if isempty(obj_idx), continue, end
    if contains(s,'pt3dadd(')
        val = str2double(regexp(s,'-?\d+(\.\d+)?','match'));
        skel.n = [skel.n; val(2:4)];
        skel.f = [skel.f; obj_idx];
        skel.d = [skel.d; val(5)];
        np = numel(skel.f); 
        if np > 1 && skel.f(end-1) == obj_idx
            skel.e = [skel.e; np-1 np];
        end

        continue
    end
    if contains(s,'}'), obj_idx = []; continue, end
    % disp(s)
  end

  if false
      %% Debug plot, check HOC correctness
      disp('debug plot of HOC cell')
      px = reshape(skel.n(skel.e,1),[],2); px(:,3) = nan;
      py = reshape(skel.n(skel.e,2),[],2); py(:,3) = nan;
      pz = reshape(skel.n(skel.e,3),[],2); pz(:,3) = nan;

      v_ = @(x) reshape(x,[],1); 

      clf, plot3(v_(px'), v_(py'), v_(pz')), axis image
      title(strrep(stub,'_','\_'))
  end

else error('Unknown file format: %s%s', stub, ext)
end

return