

function output = artificalRadon(varargin)
% radon_dat = tools.artificalRadon( RF_map, XY, ... )
% radon_dat = tools.artificalRadon( 'image', RF_map, 'xy', XY, ... )
% 
% tools.artificalRadon simulates the response of a hypothetical cell with
% a receptive field specified by RF_map and XY to a radon (flashing bar)
% protocol. 
% 
% For a usage example, look at script_EffectiveStimContrast.m
% 
% 'image', [RF_map] - receptive field of the cell to be simulated.
% 'xy', [xy_min : [] : xy_max] - spatial extent of the receptive field,
%                                in µm or degrees (coordinates of the rows
%                                and columns of [RF_map]). Assumed square.
% 'xy', [xmin xmax ymin ymax]  - alternate syntax for 'xy' (non-square).
% 
% 
% 'width',  [60]     - spacing between adjacent bars (same units as XY)
% 'angles', [6]      - 30 degree step size between orientations
% 'nBars',  [21]     - number of bars per orientation
% 'algorithm','SART' - reconstruction algorithm (see AIR-tools)
% 'stim', [struct]   - input stimulus parameters as a structure (can be
%                      more convenient than the below syntax, plus more
%                      options are exposed). 
% 
% stim_params fields:
% .step [60]     - spacing between adjacent bars (same units as XY)
% .ori_step [30] - degrees = 180/'angles' between orientations
% .n_bars [21]   - number of bars per orientation
% .width  =  1.5 * stim.step - physical width of bar on the screen
% .height = 21.5 * stim.step - physical height of bar on the screen
% 
% v0.1 - 29 August 2022 - Calvin Eiber <ceiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(n) varargin{find(named(n))+1};

%% Step 1 - set up (import/export) stimulus parameters
stim_params = any(named('stim')); 

if stim_params
  idx = find(named('stim')) + 1;
  output_params = (idx > nargin || ischar(varargin{idx}));  
  if ~output_params, stim = varargin{idx}; end
end
if ~exist('stim','var'), stim = struct; end

if any(named('width')), stim.step = get_('width');
elseif ~isfield(stim,'step'), stim.step = 60; % µm
end

if any(named('angles')), stim.ori_step = (180 / get_('angles'));
elseif ~isfield(stim,'ori_step'), stim.ori_step = 30; % degrees
end

if any(named('nbars')), nBars = get_('nbars');
elseif isfield(stim,'n_bars'), nBars = stim.nBars; 
else nBars = 21; 
end

if ~isfield(stim,'width'), stim.width = 1.5 * stim.step; end
if ~isfield(stim,'height'), stim.height = 21.5 * stim.step; end


rdat = struct; 
rdat.ori = reshape( ones(nBars,1) * (0:stim.ori_step:(180-stim.ori_step)), [], 1);
rdat.x   = reshape( linspace(-10,10,nBars)' * ones(1,numel(unique(rdat.ori))) * stim.step, [], 1); 
rdat.y   = nan(size(rdat.ori)); 

stim.ori = (-rdat.ori) * pi / 180;
stim.xpos = rdat.x .* sin(rdat.ori * pi/180);
stim.ypos = rdat.x .* cos(rdat.ori * pi/180);

if stim_params && output_params, output = stim; return, end

if any(named('image')), img = get_('image'); 
else                    img = varargin{1}; 
end

if any(named('algorithm')), rdat.algorithm = get_('algorithm');
else                        rdat.algorithm = 'sart';
end

%% Step 2 - parse model specification 

use_labels = any(named('weights')); 
use_profiles = any(named('profile'));
do_animation = any(named('animate'));

model = struct; 
if use_profiles, model.profile = get_('profile'); 
    nP = size(model.profile,2);
end
if use_labels
  idx = find(named('weights')) + 1;
  if (idx > nargin || ischar(varargin{idx})), 
    % weights specified but not followed by a structure, assume the image
    % is labelled according to "1,2,3,..." = profile_1, _2, _3, ...
    
    img(isnan(img)) = 0; 
    model.labels = unique(img(:));
    if ~isfield(model,'profile') || length(model.labels) > (1 + nP)
      warning('Treating this input as a continuous image...')
      use_labels = false; 
    else
      model.weight = zeros(length(model.labels), nP); 
      for pp = 1:nP, model.weight(model.labels == pp,pp) = 1; end
    end
    
  else
    
    cmd = varargin{idx};
    model.labels = reshape(cmd.labels,[],1);
    model.weight = cmd.weights;
  end
elseif use_profiles, nP = min(nP, size(img,3));
end

%% Get input XY span 

if any(named('xy')), roi = get_('xy'); 
else roi = varargin{2};
end

if all(size(roi) == [1 4])
    x = linspace(roi(1),roi(2),size(img,2));
    y = linspace(roi(3),roi(4),size(img,1));
elseif any(size(roi) == 2) && max(size(roi)) > 2
    if size(roi,1) < size(roi,2), roi = roi'; end
    roi = [min(roi) max(roi)];     
    x = linspace(roi(1),roi(3),size(img,2));
    y = linspace(roi(2),roi(4),size(img,1));
else
    roi = [min(roi(:)) max(roi(:))];
    x = linspace(roi(1),roi(2),size(img,2));
    y = linspace(roi(1),roi(2),size(img,1));
end
    
if use_profiles
    time = linspace(x(1),x(end),size(model.profile,1)); 
    wave = zeros(numel(time), numel(stim.ori)); 
end

[gx,gy] = meshgrid(x,y); 
d_a = 1/sqrt(numel(gy)); 

%% Generate RADON projections

if do_animation, 
    Tools.standardFigure('Name','Artificial Radon')   
    h = imagesc(x,y,img(:,:,1)); 
    axis image off xy, caxis(caxis ./ 0.8), hold on    
    if use_profiles, h(2) = plot(time, 0*time,'k-','LineWidth',2); end
end

coe = zeros(1,size(img,3));
rdat.y = rdat.y * coe;

for ii = 1:length(rdat.y)

    bar_img = abs(((gy-stim.ypos(ii)) * sin(stim.ori(ii)) + ...
                   (gx-stim.xpos(ii)) * cos(stim.ori(ii)))) ...
                          < stim.height/2 & ...
              abs(((gy-stim.ypos(ii)) * cos(stim.ori(ii)) - ...
                   (gx-stim.xpos(ii)) * sin(stim.ori(ii)))) ...
                          < stim.width/2;
    
    if use_profiles && size(img,3) > 1
      for pp = 1:nP
          i_plane = img(:,:,pp);
          coe = sum(i_plane(bar_img)) * d_a;
          wave(:,ii) = wave(:,ii) + model.profile(:,pp) * coe;
      end
    elseif use_labels
      coe = arrayfun(@(n) sum(img(bar_img) == n), model.labels) * d_a; 
      coe = (model.weight' * coe);
      wave(:,ii) =  model.profile * coe;
    else
      for pp = 1:size(img,3)
        coe(pp) = sum(img(bar_img)) * d_a;
        rdat.y(ii,:) = coe;
      end
    end
   
    %% Update figure if relevent
    if ~do_animation || ~any(ishandle(h)), continue, end    
    h(1).CData = (0.2*bar_img + img(:,:,1));
    if use_profiles, h(2).YData = wave(:,ii);
    else title(sprintf('%0.5f',coe(1)))
    end
    pause(0.1)
end

%% Compute inverse radon transform of input model 

if ~exist('wave','var') 
    rdat.y_all = rdat.y;
    if any(named('no-fig')), output = rdat; return, end
    
    plots.standardFigure('Name','Artificial Radon')
    rdat = plots.plot_radon_IMG(rdat);
    output = rdat; return
end

%% Compute PCA of simulated temporal profiles 

wave = wave / max(wave(:));

if any(named('noise')), noise = get_('noise');
  if ~all(size(noise) == size(wave))
      noise = randn(size(wave)) * noise(1); 
  end
else noise = 0;
end

[coeff,score] = pca(wave + noise,'Centered',false);

if any(named('nK')), nK = get_('nK');
else                 nK = 4; 
end

coeff = coeff(:, 1:nK);
score = score(:, 1:nK); 
coe_scale = max(coeff);

coeff = coeff * diag(1./coe_scale);
score = score * diag(coe_scale);

rdat.y_all = coeff; 
rdat.wave = score;

rdat = Tools.plot_radon_IMG(rdat);

h = flipud(get(gcf,'Children')); 

for kk = 1:nK
        
    this = h(3*kk-2); 
    x = linspace(this.XLim(1),this.XLim(end),numel(time));    
    y = (0.8*score(:,kk)+1)*range(this.YLim)/2;
    
    hold(this,'on')
    plot(x,y,'w-','LineWidth',1.5,'Parent',this)
    
end

output = rdat; 

end

