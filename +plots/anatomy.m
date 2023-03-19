function anatomy(anat, varargin) 
% plots.anatomy(anat_data, [radon_dat], ...)
% plots.anatomy(anat_data, [radon_dat], [gauss_model], ...)
% 
% Generate simple heatmap showing 2D alignment of cell and RF. 
% If [radon_dat] is not supplied, only the anatomy is plotted and most
% options are disabled. 
% 
% Options: 
% -xy [x y]  : set cell center in radon heatmap. If not set, the user is
%               prompted to select the XY center. 
% -gm [gm]   : set gaussian model (optional, overrides positional argument)
% -clf       : suppress figure clearing (enables to write on existing axes)
% -zoom [10] : set zoom-in factor on cell in radon image
% -date []   : set experiment date (for alignment of cell and heatmap),
%              defaults to the date in anat.name (e.g. 20210526)
% -id [1]    : select image ID to view
% -g_id [1]  : select gaussian contour to plot (if relevent)
% -recenter  : move fitted gaussian centers to selected cell position
% 
% v0.2 - 23 December 2023 - Calvin Eiber <c.eiber@ieee.org>
% adapted from v1 code (Elissa Belluccini)

if nargin > 1 && isstruct(varargin{1})
     rdat = varargin{1};
else rdat = []; 
end

if nargin > 2 && isstruct(varargin{2}) && ~ischar(varargin{1}) 
     gm = varargin{2};
else gm = []; 
end

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

% Step 1, get anatomy 2D contours 

anat_x = reshape(anat.node(anat.edge,1),[],2); 
anat_x = anat_x - median(anat.node(:,1)); 
anat_x(:,3) = nan;

anat_y = reshape(anat.node(anat.edge,2),[],2); 
anat_y = anat_y - median(anat.node(:,2));
anat_y(:,3) = nan;

% anat_z = reshape(anat.node(anat.edge,3),[],2); anat_z(:,3) = nan;

v_ = @(x) reshape(x,[],1); 

if ~any(named('-clf')), clf, end

if ~isempty(rdat)

  if any(named('-zoom')), subplot(1,2,1); end
  
  img_id = 1; 
  if any(named('-id')), img_id = get_('-id'); end
  img = rdat.images{img_id};
    
  % get experimental date from anat.filename
  exp_date = regexp(anat.name,'20\d{6}','match','once');
  if any(named('-date')), exp_date = get_('-date'); end
  if ischar(exp_date), exp_date = str2double(exp_date); end
  if isnan(exp_date) || isempty(exp_date), exp_date = 99999999; 
      warning('Experiment date in "%s" invalid or not found', anat.name);
  end

  % For hoc files in the wrong orienation: Orient RF map to match cell
  % orientation (in conjection with correction from 'Orienting Cell Morphology and RF Map.docx')
  if strcmp(anat.name,'20190904_Cell_02')
      h_im = imagesc(rdat.range,rdat.range,imrotate(img,101.17));
  elseif strcmp(anat.name,'20200116_Cell_02')
      h_im = imagesc(rdat.range,rdat.range,imrotate(img,187));
  elseif strcmp(anat.name,'20210601_Cell_01')
      h_im = imagesc(rdat.range,rdat.range,flipud(imrotate(img,-250)));
  elseif any(strcmp(anat.name,{'20211122_Cell_02','20211122_Cell_02','20211122_Cell_04'}))
      h_im = imagesc(rdat.range,rdat.range,flipud(imrotate(img,-265.25)));

  % See PRJ-VisionLab\EAB\Notes\Orienting Cell Morphology and RF Map.docx
  % ASB correction: Rotate image 90 deg anticlockwise and flip about 
  %                  horizontal AND vertical axis
  elseif (exp_date < 20210526)
      % TODO - can we please replace the call to imrotate with a call to
      % rot90 (much faster)? 
    h_im = imagesc(rdat.range,rdat.range,rot90(img));
  else
  % MFB correction: Rotate image 90 deg anticlockwise and flip about 
  %                  horizontal axis only
    h_im = imagesc(rdat.range,rdat.range,flipud(rot90(img,-1)));       
  end

   % Crop extra pixels after rotation:
  if length(h_im.CData)>101
     excess = (length(h_im.CData)-101);
     trim1 = floor(excess/2);
     trim2 = ceil(excess/2);
     h_im.CData = h_im.CData(trim1+1:end-trim2,trim1+1:end-trim2);
  end
  axis square, axis image off xy, hold on

  if any(named('-redblue'))
    set(gca,'CLim',[-1 1] * max(abs(get(gca,'CLim'))))
    rbcm = interp1((-5:5)', redbluecmap, linspace(-5,5,101)); 
    colormap(gca,rbcm)
  end

  if any(named('-xy')), mid_xy = get_('-xy');
  else
    disp('Select receptive field centre');
    [mid_xy(1),mid_xy(2)] = ginput(1);
    fprintf('-XY: [%0.2f %0.2f]\n', mid_xy)
  end

  plot(v_(anat_x' + mid_xy(1)), ...
       v_(anat_y' + mid_xy(2)), 'k-','LineWidth',1.3)

  cb = colorbar; 
  cb.TickDirection = 'out';
  cb.Box = 'off'; 
  cb.Location = 'westoutside';
  cb.Position = [0.055,0.11,0.025,0.8]; % hard-coded 

  axis image on xy
  try tidyPlotForIllustrator, end
  
  if any(named('-zoom'))

    try zf = get_('-zoom'); % percentage padding
    catch err, zf = 10; %#ok<NASGU> 
    end
    if ischar(zf), zf = 10; end

    zf = zf / 100; 

    anat_lim = [min(anat_x(:), [], 'omitnan') ...
                max(anat_x(:), [], 'omitnan') ...
                min(anat_y(:), [], 'omitnan') ...
                max(anat_y(:), [], 'omitnan') ];
      
    new_xy = anat_lim * [1+zf -zf 0 0; 
                         -zf 1+zf 0 0; 
                         0 0 1+zf -zf; 
                         0 0 -zf 1+zf];

    [~,i1] = min(abs(rdat.range-new_xy(1))); %i1 = i1-1;
    [~,i2] = min(abs(rdat.range-new_xy(2))); %i2 = i2+1;
    [~,i3] = min(abs(rdat.range-new_xy(3))); %i3 = i3+1;
    [~,i4] = min(abs(rdat.range-new_xy(4))); %i4 = i4-1;

    subplot(1,2,2);
    imagesc(rdat.range(i1:i2),rdat.range(i3:i4),h_im.CData(i3:i4,i1:i2));
    axis image on xy, hold on
    
    plot(v_(anat_x' + mid_xy(1)), ...
         v_(anat_y' + mid_xy(2)), 'k-','LineWidth',1.3)

    try tidyPlotForIllustrator, end %#ok<TRYNC> 
  end

else
  % just show the cell anatomy

  mid_xy = [0 0];
  if any(named('-xy')), mid_xy = get_('-xy'); end
  plot(v_(anat_x' + mid_xy(1)), ...
       v_(anat_y' + mid_xy(2)), 'k-','LineWidth',1.3)
  hold on
  try tidyPlotForIllustrator, end %#ok<TRYNC> 

end

%% Add Gaussian radii to model 

% TODO - 


if any(named('-gm')) || ~isempty(gm)

  if any(named('-gm')), gm = get_('-gm'); end

  C = lines(gm.n_gaussians);
  style = {'LineWidth',1.2,'Color'};

  for gg = 1:gm.n_gaussians

    theta = linspace(0,2*pi,91);
    oneSD_circle = [cos(theta); sin(theta)]'; 

    if isfield(gm,'gauss_eccentricity')
      ecc = sqrt(1 - (gm.gauss_eccentricity(gg) .* ...
                  cos(theta-deg2rad(gm.gauss_angle(gg))).^2));
      oneSD_circle = [cos(theta)./ecc; sin(theta)./ecc]';
    end

    xy = gm.center_xy(gg,:) + gm.gauss_radius(gg) * oneSD_circle;

    % THESE XY values are in the original (not rotated or flipped)
    % coordinate frame. this needs correction 

    % See PRJ-VisionLab\EAB\Notes\Orienting Cell Morphology and RF Map.docx
    % ASB correction: Rotate image 90 deg anticlockwise and flip about 
    %                  horizontal AND vertical axis
    if (exp_date < 20210526)
          % TODO - can we please replace the call to imrotate with a call to
          % rot90 (much faster)? 
      % h_im = imagesc(rdat.range,rdat.range,rot90(imrotate(img,-90),2));
      xy = xy(:,[2 1]) .* [1 -1];
    else
    % MFB correction: Rotate image 90 deg anticlockwise and flip about 
    %                  horizontal axis only
      xy = xy(:,[2 1]) .* [1 1];
      warning('TODO please check correct rotation')
      % h_im = imagesc(rdat.range,rdat.range,flipud(imrotate(img,-90)));       
    end


    if any(named('-recenter'))
        xy = xy + mid_xy - gm.center_xy(1,:);
    end

    plot(xy(:,1), xy(:,2), style{:}, C(gg,:))
  end
end

