
function f = standardFigure(varargin)
% Creates a figure of a standard size, and optionally loads a couple of
% utility tools into the calling workspace. 
% 
%   EXAMPLE USAGE: 
% Tools.standardFigure('1column','Name','my-fig-name','Height',1.2,'Tools')
%   If a figure is named "my-fig-name", use it and clear it. Otherwise,
% make a new figure with that name sized for 1 column with 1.2x height.
% Additionally, return the standard tools (G,C,W). 
%
% Named figure sizes: [Portrait, Landscape, 1column, 2column, 1.5column]

named = @(n) strncmpi(varargin,n,3);

% If requested, add tools to calling workspace
if any(named('Tools'))    
    assignin('caller','G',@(v) [v v v]/10); 
    assignin('caller','C',lines(7)); 
    evalin('caller','W = @(i,v)(C(mod(i-1,7)+1,:)+max(0,v))/(1+v);')
    evalin('caller','tweak = @(h,p) set(h,''Position'',h.Position + p/100);');
end

% Use a persistent named singleton figure
if any(named('Name'))
    figure_name = varargin{find(named('Name'),1)+1};
    f = findobj('Name',figure_name);
    if ~isempty(f)        
        figure(f); clf 
        if nargout == 0, clear, end
        return
    end    
end

f = figure('Color','w');
if any(named('Name')), f.Name = figure_name; end

dP = [0 0 0 0];

% Preset values for position
if any(named('Portrait')),      f.Position = [10 10 820 980];
elseif any(named('Landscape')), f.Position = [40 100 1200 820];
elseif any(named('1column')),   dP = [0 0 0 0];
elseif any(named('2column')),   dP = [-0.3 0 0.5 0];
elseif any(named('1.5column')), dP = [-0.1 0 0.25 0];
end

if any(named('Height')), 
    dP([2 4]) = [-1 1] * (varargin{find(named('Height'),1)+1} - 1);
end

f.Position = f.Position + f.Position([3 4 3 4]).*dP;

if nargout == 0, clear, end


end