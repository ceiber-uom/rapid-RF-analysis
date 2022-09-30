
function data = inverseRadon(data)
% 
% 
% 
% 
% Revised iRadon 21 September 2018 uses AIR-Tools-II to compute inverse
% radon results... faster and more accurate than rolling my own. CDE

if ~isfield(data,'algorithm'), data.algorithm = @sart; end
if ischar(data.algorithm), data.algorithm = str2func(data.algorithm); end
if isempty(which(func2str(data.algorithm))),    
%     path(path,locateFiles(['\\research-data.shared.sydney.edu.au\SMS\' ... 
%                                    'PRJ-VNRG\PHYSIOL\MFILES\AIR-TOOLS']));
    path(path,'../HEKA Radon/AIR-Tools');
    run('AIRToolsII_setup.m')
end

if ~isfield(data,'ori'), data = utils.prepareRadon(data); end








data.range  = linspace(min(data.x),max(data.x),101); % <<<<< CHANGE RESOLUTION HERE

A = calc_system_matrix(data); % Sparse system matrix

data.system_matrix = A; 

if isfield(data,'AIR_opts'), opts = data.AIR_opts; 
else                         opts = struct; 
end

nX = numel(data.range); 
nY = numel(data.ori); 
nA = numel(unique(data.ori));

if isfield(data,'iterations'), nIter = data.iterations; 
else                           nIter = 50; 
        opts.stoprule.type = 'NCP';
        opts.stoprule.res_dims = [nY/nA nA];
end

if ~isfield(data,'censor') || data.censor
     Y = data.y .* sqrt(full(sum(A,2) / max(sum(A,2)))); % Censor partially hidden bars
else Y = data.y;
end

if strcmp(func2str(data.algorithm),'fbp')
    
    Y = Y / mean(A(:)) / sqrt(2); 
    
    angles = deg2rad(unique(data.ori)); 
    if isfield(data,'filter')
         image = data.algorithm(A,Y,angles,data.filter);
    else image = data.algorithm(A,Y,angles);
    end
else    
    image = data.algorithm(A,Y,nIter,[],opts);
end

% Meta-data about this fit
data.image = reshape(image, nX, nX); 

[gX,gY]  = meshgrid(data.range,data.range);

% Include only valid centre region
inc = sqrt(gX.^2 + gY.^2);
inc = inc <= gX(end); 
data.mask = inc; 

if ~isfield(data,'fit'), return, end
%%

% % identify receptive field
% field = reshape(field, size(gX));
% 
% if all(field(:) == 0), return, end
% 
% if options.doComplexRadon
%     [~,idx] = max(abs(field(:).*inc(:)));
%     % Normalize so maximum falls at [0i + 1*(scale)]. 
%     field = real(field ./ field(idx) * abs(field(idx)));    
%     % If we were to be naive and just take abs(field), the oscillations
%     % would make it too hard to fit. 
% end
% 
% gauX  = '( cos(deg2rad(phi)) * (x-xc) - sin(deg2rad(phi)) * (y-yc) )';
% gauY  = '( sin(deg2rad(phi)) * (x-xc) + cos(deg2rad(phi)) * (y-yc) )';
% 
% gau2d = sprintf('a * exp( - (%s/sx)^2 - (%s/sy)^2 ) + c', gauX, gauY);
% 
% ftype            = fittype( gau2d, ...
%                    'independent', {'x', 'y'}, 'dependent', 'z' );
% fopts            = fitoptions( ftype );
% fopts.Display    = 'Off';
% 
% [maxR, maxI]     = max(field(:).*inc(:));
% [minR, ~   ]     = min(field(inc));
% [maxY, maxX]     = ind2sub(size(field), maxI);
% vfx_max          = vRange(maxX);
% vfy_max          = vRange(maxY);
% 
% maxSD = max(gX(:)); 
% 
% 
% minVF            = vRange(1);
% maxVF            = vRange(end);
% fopts.Lower      = [ -Inf  -Inf  -90  0       0        minVF    minVF ];
% fopts.StartPoint = [ maxR  minR    0  maxSD/3 maxSD/5  vfx_max  vfy_max ];
% fopts.Upper      = [  Inf   Inf   90  maxSD   maxSD    maxVF    maxVF ];
% 
% [rfFit, rfGoF] = fit( [gX(inc), gY(inc)], field(inc), ftype, fopts );
% 
% % Place output in structure
% data.RF_fit = rfFit;
% data.RF_stats  = rfGoF;

end

function field = OLD_radon(varargin) %#ok<DEFNU>

data = evalin('caller','data'); % Cheat, we're basically ignoring everything

% image = data.algorithm(A,Y,nIter,[],opts);

angles = unique(data.ori);

% Get mean respone for each cell at each ori
X = arrayfun(@(ori) data.x(data.ori == ori), angles, 'UniformOutput', false);
Y = arrayfun(@(ori) data.y(data.ori == ori), angles, 'UniformOutput', false);
for ii = 1:length(Y), Y{ii}(isnan(Y{ii})) = 0; end

cellget = @(varargin) cellfun(varargin{:}, 'Unif',0);


%   used hamming window 0.6 of sampling frequency
% uncomment below to use equivalent butterworth, easier to program
% need to evaluate most appropriate filter with real noisy data
[filtB,filtA] = butter(1,0.8);
Y_filt = cellget(@(y) filtfilt(filtB, filtA, y), Y);

% inverse radon transform - backprojection method
 vRange  = linspace(min(data.x),max(data.x),101);
[gX,gY]  = meshgrid(vRange,vRange);

angles = deg2rad(angles);
field = nan(numel(gX), numel(angles));

for ii = 1:numel(angles),
    tX = gX(:) * sin(angles(ii));
    tY = gY(:) * cos(angles(ii));
    field(:,ii) = interp1(X{ii}(:),  Y_filt{ii}(:), tX+tY, 'pchip',NaN);
end

field = nanmean(field,2); 

end

function A = calc_system_matrix(data)


if isfield(data,'info'), info = data.info; else info = struct; end

angles = unique(data.ori);
X = arrayfun(@(ori) data.x(data.ori == ori), angles, 'UniformOutput', false);    
dX = median(cellfun(@(x) mean(diff(x)), X)); 

nX = numel(data.range); 
nY = numel(data.ori); 
nA = numel(angles); 
nB = numel(unique(data.x));

[gY,gX]  = meshgrid(data.range,data.range);

if ~isfield(info,'width'), info.width = 1.5; end % warning('Usually 1.5')
if ~isfield(info,'height'), info.height = (nB-1) + info.width; end

stim.width = info.width * dX;
stim.height = info.height * dX;

% CHECK these values against expoDataSet values!

stim.ori = deg2rad(data.ori + 90); 
stim.xpos =  sin(stim.ori) .* data.x;
stim.ypos = -cos(stim.ori) .* data.x;

bar_area = (stim.width .* stim.height); % units of um²
px_area = mean(diff(data.range)) .^ 2; % units of µm²

persistent system_matrix system_info

use_cached = ~isempty(system_info) && ~isempty(system_matrix) && ...
                  all(system_info == [nX nY nA info.width info.height]); 

if use_cached
    % disp('Using cached system matrix ... ')
    A = system_matrix * (px_area / bar_area);
    return
end

system_info = [nX nY nA info.width info.height]; 

%% Build system matrix (taking into account width of stimuli)

% fundamental equation: 
% A * IMG = DATA      (A = 126 x 10201 system matrix)
% We could directly compute this, but it's sensitive to noise (ill-posed)

% (is it worth adding a wee bit of gaussian blur here?)

idx = [];
jdx = [];
pix = [];

[ssY,ssX] = meshgrid(linspace(-1,1,5)/2 * mean(diff(gX(:,1))));
radius = inf; % max(gX(:))^2;

printInfo;

for ii = 1:numel(gX)
    
    if mod(ii,21) || ii == numel(gX)
        printInfo('Building System Matrix (%0.1f%%) ', 100*ii/numel(gX))
    end
    
    if (gX(ii)^2 + gY(ii)^2) > radius, continue, end
    
    row = false(nY,numel(ssX));
    
    sin_ori = sin(stim.ori);
    cos_ori = cos(stim.ori);
    
    pix_y = gY(ii) - stim.ypos;
    pix_x = gX(ii) - stim.xpos;
    
    for ss = 1:numel(ssX)
        row(:,ss) = abs((pix_y+ssY(ss)) .* sin_ori + ...
                        (pix_x+ssX(ss)) .* cos_ori) < stim.height/2 & ...
                    abs((pix_y+ssY(ss)) .* cos_ori - ...
                        (pix_x+ssX(ss)) .* sin_ori) < stim.width/2 ;
    end
    
    row = sum(row,2)/ss; 
    rix = find(row); 
    
    % ensure that A contains a value for each angle.
    if numel(unique(stim.ori(rix))) < nA, continue, end
    
    idx = [idx; ii+0*rix];  %#ok<AGROW>
    jdx = [jdx; rix];       %#ok<AGROW>
    pix = [pix; row(rix)];  %#ok<AGROW>
end

A = sparse(jdx,idx,pix,nY,numel(gX));
system_matrix = A; 

A = A * (px_area / bar_area); 
disp('Done!')

% %%
% return
% figure %#ok<UNRCH>
% h = imagesc(gX(:,1), gX(:,1), full(reshape(A(1,:),nX,nX)));
% axis image off xy
% 
% for ii = 1:nY    
%     h.CData(:) = A(ii,:);
%     pause(0.1)
% end

end


