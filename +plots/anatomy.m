function anatomy(anat, varargin) 
% plots.anatomy(anat_data, [radon_dat], ...)

if nargin > 2 && isstruct(varargin{1})
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

anat_x = reshape(anat.node(anat.edge,1),[],2); anat_x(:,3) = nan;
anat_y = reshape(anat.node(anat.edge,2),[],2); anat_y(:,3) = nan;
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

  % See notes: PRJ-VisionLab\Elissa_Belluccini\Notes\
  %                             Orienting Cell Morphology and RF Map.docx
  % ASB correction: Rotate image 90 deg anticlockwise and flip about 
  %                  horizontal AND vertical axis
  if (exp_date < 20210526)
    h_im = imagesc(rdat.range,rdat.range,rot90(imrotate(img,-90),2));
  else
  % MFB correction: Rotate image 90 deg anticlockwise and flip about 
  %                  horizontal axis only
    h_im = imagesc(rdat.range,rdat.range,flipud(imrotate(img,-90)));       
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
  end

  plot(v_(anat_x + mid_xy(1)), ...
       v_(anat_y + mid_xy(2)), 'k-','LineWidth',1.3)

  cb = colorbar; 
  cb.TickDirection = 'out';
  cb.Box = 'off'; 
  cb.Location = 'westoutside';
  cb.Position = cb.Position .* [0.3,0.75,1.2,1.82]; % hard-coded 

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
    imagesc(rdat.range(i1:i2),rdat.range(i4:i3),h_im.CData(i4:i3,i1:i2));
    axis image on xy, hold on
    
    plot(v_(anat_x + mid_xy(1)), ...
         v_(anat_y + mid_xy(2)), 'k-','LineWidth',1.3)

    try tidyPlotForIllustrator, end
  end

else
  % just show the cell anatomy

  mid_xy = [0 0];
  if any(named('-xy')), mid_xy = get_('-xy'); end
  plot(v_(anat_x + mid_xy(1)), ...
       v_(anat_y + mid_xy(2)), 'k-','LineWidth',1.3)
  hold on
  try tidyPlotForIllustrator, end

end



if any(named('-gm')) || ~isempty(gm)

    g_id = 1; 
    if any(named('-gm')), gm = get_('-gm'); end
    if any(named('-gid')), g_id = get_('-gid'); end

    C = lines(max(7,g_id));

    theta = linspace(0,2*pi,91);
    oneSD_circle = [cos(theta); sin(theta)]'; 

    if isfield(gm,'gauss_eccentricity')
      ecc = sqrt(1 - (gm.gauss_eccentricity(gg) .* ...
                  cos(theta-deg2rad(gm.gauss_angle(gg))).^2));
      oneSD_circle = [cos(theta)./ecc; sin(theta)./ecc]';
    end

    style = {'-','Color',C(g_id,:),'LineWidth',1.2};
    % xy = gm.center_xy(g_id,:) + gm.gauss_radius(g_id) * oneSD_circle;
    xy = mid_xy + gm.gauss_radius(g_id) * oneSD_circle;
    plot(xy(:,1), xy(:,2), style{:})

end

