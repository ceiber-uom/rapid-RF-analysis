function anat = loadAnatomy( filename, varargin )
% anat = loadAnatomy( filename, ... )
% 
% Load a tracing of a cell from one of the supported file formats and
% return a structure containing the dendritic tree anatomy for that cell.
% If no return argument is specified, a structure called 'anat' is
% automagically created in the calling workspace. 
% 
% This code is comfortable handling three kinds of input data: 
% - .mat files corresponding to the skeletonised retinal ganglion cells 
%        published in Bae et al. (2018)
%       (https://doi.org/10.1016/j.cell.2018.04.040)
%       (https://github.com/seung-lab/e2198-gc-analysis) 
% - .mat files traced in matlab, published by Eiber et al. (2019)
% - .hoc files exported from traced RGCs using Imaris 
% 
% Returned dendrites are in units of µm (nominally) and are returned as:
% anat.nodes : list of nodes (points in xyz) 
% anat.edges : list of edges between nodes (index into anat.nodes)
% 
% if available, the dendrite diameter (in µm) and soma location are also
%               returned.
% 
% Issues around tissue shrinkage, realignment, etc. are not handelled.
% 
% Options:
% -scale [s]        : apply scale factor (Default: 1 for .hoc and Eiber
%                      (2019) files, [0.092 0.066 0.066] for Bae (2018)). 
% -scale [sx sy sz] : apply non-uniform scale (e.g. anisotropic voxels)
% -single-tree : unify results into a single dendritic tree (i.e. connect
%                primary dendrites). Requires knowledge of soma location.
% -plot : show the loaded anatomy in a 3D plot. Can also be used to show 
%         again loaded data, e.g. tools.loadAnatomy(anat, '-plot')
% 
% v0.1 - 17 September 2022, Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

if nargin < 1, filename = '?'; end
if isstruct(filename)
  if any(named('-plot'))
    anat = filename; validation_plot(anat), return, end
end
if ~exist(filename,'file'), filename = pick_file; end

do_plot = any(named('-plot'));
anat = struct;

% Determine how to correctly parse file based on filename
[~,stub,ext] = fileparts(filename);

if strcmp(ext,'.mat')

    vars = whos('-file',filename);
    
    if any(strcmp({vars.name},'dendrites'))

        disp(['Loading ' stub ext])
        in = load(filename,'dendrites','soma');

        %% convert NAN-seperated vector to nodes+edges representation
        [xy,row_id,node_id] = unique(in.dendrites,'rows');
        fil_id = cumsum(isnan(in.dendrites(:,1))); 
        
        row_id(isnan(xy(:,1))) = []; 
        xy(isnan(xy(:,1)),:) = []; 
                
        anat = struct; 
        anat.name = stub;
        anat.node = [xy 0*xy(:,1)]; 
        anat.edge = []; 

        for ii = 2:numel(node_id)
            if any(isnan(in.dendrites(ii-[0 1]))), continue, end            
            anat.edge = [anat.edge; reshape(node_id(ii-[0 1]),1,2) ];
        end

        anat.soma = [median(in.soma) 0];
        anat.f_id = fil_id(row_id); 
        anat.scale = 1;

        if any(named('-single-tree'))
            %% Add connection from each primary dendrite to nominal soma
            s_radius = sqrt(mean(sum((in.soma - anat.soma(:,1:2)).^2,2)));

            G = graph(anat.edge(:,1), anat.edge(:,2));
            pid = G.conncomp;
            s_node = size(anat.node,1) + 1;

            for pp = 1:max(pid)
                sel = find(pid == pp);
                [dist,idx] = min(sum((anat.node(sel,:) - anat.soma).^2,2));
                if sqrt(dist) > 2*s_radius, continue, end
                anat.edge = [anat.edge; sel(idx) s_node];
            end
            anat.node(end+1,:) = anat.soma;
        end

        clear in row_id fil_id node_id xy G pid s_node

    elseif all(ismember({vars.name},'root','n','e'))
        disp(['Loading ' stub ext])
        in = load(filename,'n','e'); 

        anat.name = stub;
        anat.node = in.n;
        anat.edge = in.e;
        % todo ... determine f_id, diam? 

        anat.scale = [92 66 66] / 1e3; % from the readme
        % https://github.com/seung-lab/e2198-gc-analysis

    else error('unknown internal format: %s.mat', stub)
    end

elseif strcmp(ext,'.hoc') || strcmp(ext,'.geo')
  %% Parse HOC file
  disp(['Reading ' stub ext])
  fi = fopen(filename,'rt');
  oc = onCleanup(@() fclose(fi));
  
  anat.name = stub;
  anat.node = []; % node XYZ
  anat.edge = []; % skeleton edges in graph
  anat.f_id = []; % filament ID
  anat.diam = []; % filament diameter (not used)
  anat.scale = 1;

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

        if val(1), val(1) = find(anat.f_id == oid(1),1,'last');
        else       val(1) = find(anat.f_id == oid(1),1);
        end
        if val(2), val(2) = find(anat.f_id == oid(2),1,'last');
        else       val(2) = find(anat.f_id == oid(2),1);
        end
        anat.edge = [anat.edge; val];
        continue
    end

    if isempty(obj_idx), continue, end
    if contains(s,'pt3dadd(')
        val = str2double(regexp(s,'-?\d+(\.\d+)?','match'));
        anat.node = [anat.node; val(2:4)];
        anat.f_id = [anat.f_id; obj_idx];
        anat.diam = [anat.diam; val(5)];
        np = numel(anat.f_id); 
        if np > 1 && anat.f_id(end-1) == obj_idx
            anat.edge = [anat.edge; np-1 np];
        end

        continue
    end
    if contains(s,'}'), obj_idx = []; continue, end
    % disp(s)
  end

else error('Unknown file format: %s%s', stub, ext)
end


if any(named('-sc')),  anat.scale = get_('-sc');
elseif any(named('-vo')), anat.scale = get_('-vo')/1e3;
end

anat.node = anat.node .* anat.scale; 
if do_plot || nargout == 0, validation_plot(anat), end
if nargout == 0, assignin('caller','anat',anat); clear, end

return


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



function validation_plot(anat)
%% Debug plot, check loaded anatomy correctness

px = reshape(anat.node(anat.edge,1),[],2); px(:,3) = nan;
py = reshape(anat.node(anat.edge,2),[],2); py(:,3) = nan;
pz = reshape(anat.node(anat.edge,3),[],2); pz(:,3) = nan;

v_ = @(x) reshape(x,[],1); 

cla, plot3(v_(px'), v_(py'), v_(pz')), axis image, hold on

if isfield(anat,'soma'), 
  plot3(anat.soma(:,1),anat.soma(:,2),anat.soma(:,3),'o')
end
if isfield(anat,'name'), title(strrep(anat.name,'_','\_')), end
