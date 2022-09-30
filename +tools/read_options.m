
function options = read_options( arg_in, varargin )
% function options = setup_options ( input arguments, default options )
% 
% This utility function accepts the 'varargin' cell array from the calling
% function and the list of default options, and returns a structure array
% of options. 
% 
% Options may be specified by passing in an options array (optionally
% preceeded with the -opt flag) or by key-value pairs. 


if nargin == 0, error('not enough input arguments'), end

if nargin > 1 && iscell(varargin{1}), defaults = varargin{1}; 
else defaults = varargin; 
end

if mod(numel(defaults),2), error('please use key-value pairs'); end
if all(cellfun(@ischar, defaults(:,1))) && size(defaults,2) == 2, 
    defaults = defaults'; % convert from tall array
end

options = struct; 

% note arg_in is the vararin{} of the calling workspace 
named = @(v) strncmpi(v,arg_in,length(v));
get_ = @(v) arg_in{find(named(v))+1};


find_options = cellfun(@isstruct,varargin);        
if any(named('-opt')), options = get_('-opt');
elseif any(find_options), options = varargin{find(find_options,1)}; 
end

abbr_keys = cell(size(defaults)); 

for ii = 1:2:numel(defaults)

    for ss = 1:numel(defaults{ii})
      scv = cellfun(@(d) strncmpi(defaults{ii},d,ss), defaults(1:2:end));
      if sum(scv) == 1, break, end
    end

    abbr_keys{ii} = ['-' lower(defaults{ii}(1:ss))];

    if ~isfield(options, defaults{ii}),
        options.(defaults{ii}) = defaults{ii+1};
    end
    
    if any(named(defaults{ii})) % was it specified by name? 
        options.(defaults{ii}) = get_(defaults{ii});
    elseif any(named(abbr_keys{ii})) % or by shorthand?
        options.(defaults{ii}) = get_(abbr_keys{ii});
    end

    % interpret as [partial] path, these need to end in [filesep] 
    if ischar(defaults{ii+1}) && ~isempty(defaults{ii+1}) && ...
                                 ismember(defaults{ii+1}(end),'/\') && ...
                                ~isempty(options.(defaults{ii})) && ... 
                                ismember(options.(defaults{ii})(end),'\/_')
      options.(defaults{ii})(end+1) = filesep; 
    end
end

if any(named('--list-options'))
  for ii = 1:2:numel(defaults)
    fprintf('%s (%s): ', defaults{ii}, abbr_keys{ii})
    disp(options.(defaults{ii}))
  end
  error('options list displayed.')
end