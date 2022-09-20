
function density_image = dendriticDensity( anat, varargin )
% img = dendriticDensity( anatomy, [radon data], ... )
% 
% Create a 2D dendritic density image of the supplied cell anatomy in the
% coordinate system of the supplied radon data. works with anat objects
% returned by tools.loadAnatomy. If no input arguments are supplied, the
% user is prompted to select a cell. 
% 
% if [radon data] is not supplied, some default coordinate system based on
% the extent of the cell anatomy is used. The resolution of this default
% coordinate system can be controlled using the -n or -resol switches. 
%
% IT IS ASSUMED THAT THE RADON COORDINATE SYSTEM AND THE ANATOMY COORDINATE
% SYSTEM USE THE SAME UNITS (MICROMETERS). IF THIS IS NOT THE CASE, THEY
% MUST BE CONVERTED TO THE SAME UNITS. 
% 
% Options:
% -zero [x y] : the coordinates of (0,0) in the anatomy correspond to what 
%               coordinates in the RDAT coordinate system?
% -center     : center the cell in the RDAT coordinate system
% -plot       : make visualisation of 
% 
% -align      : TO BE IMPLEMENTED...
% 
% -resol [x]  : set grid resolution to x (presumably um).
%               Has no effect if radon data supplied.
% -n [50]     : set grid resolution to n (fraction of range of X).
%               Has no effect if radon data supplied. 
% -x [xmax]   : set data extent to x (default: exactly conver the cell). 
%               Has no effect if radon data supplied.
% 
% v2.0 - 19 September 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

if nargin == 0, anat = tools.loadAnatomy; end

if any(named('-im')), rdat = get_('-im'); % e.g. -image [radon data]
elseif nargin > 1 && isstruct(varargin{1}), rdat = varargin{1};
else rdat = []; 
end

if any(named('-z')), xy_zero = get('-z');
elseif any(named('-cent')), xy_zero = -median(anat.node(:,1:2));
else xy_zero = [0 0];
end

if isempty(rdat)

    xmax = max(max(abs(anat.node(:,1:2) + xy_zero))); 
    if any(named('-x')), xmax = get_('-x'); end

    do_warning = false;
    if any(named('-re')), xres = get_('-re'); 
    elseif any(named('-n')), xres = xmax / get_('-n');
    else xres = xmax/50; do_warning = true;
    end
    
    rdat.range = 0:xres:xmax;
    rdat.range = reshape(unique([-rdat.range rdat.range]),1,[]);
    
    if do_warning
        warning('no RDAT supplied, using [%dx%d] image, +-%g', ...
                 numel(rdat.range),numel(rdat.range),xmax)
    end
end

mk_image([]); % reset

if any(named('-align'))

    error TODO_alignment_code_from_choice_15
end

density_image = mk_image(rdat, anat, xy_zero);

if nargout == 0 || any(named('-p'))
  show_density(rdat, anat, density_image, xy_zero)
end

if nargout == 0, clear, end
return


%% Compute the image from the tracing
function img = mk_image(rdat,anat,xy_zero)

if nargin < 3, xy_zero = [0 0]; end
if nargin < 2, anat = evalin('caller','anat'); end
if nargin < 1, rdat = evalin('caller','rdat'); end

% this may be called in an optimisation loop ... precompute variables
persistent units
if isempty(rdat), units = []; return, end % reset
if isempty(units)
    units = (1:length(rdat.range))  / [ rdat.range; 0*rdat.range+1 ];
end

%%
dend_xy = (anat.node(:,1:2) + xy_zero)*units(1) + units(2) + 0.5;
% pixelized co-ordinates: a point is in px(x=1) if 1 < p.x < 2

img = zeros(numel(rdat.range), numel(rdat.range)); 

for dd = 1:size(anat.edge,1)    
    
    xy = dend_xy(anat.edge(dd,:),:);
    
    while true % propegate along dendrite adding to each touched voxel
        
        dy = diff(xy);
        ij = floor(xy(1,:));        
        
        if any(ij < 1 | ij([2 1]) > size(img)), break, end
        
        r = (ij([1 1 2 2])-xy([1 1 3 3])+[0 1 0 1])./dy([1 1 2 2]); 
        v = min(r(r>0)); 
        
        if v > 1
            img(ij(2),ij(1)) = img(ij(2),ij(1)) + sqrt(sum(dy.^2)); 
            break
        elseif r(1) == v, mxy = [ij(1)    r(1)*dy(2)+xy(3)];
        elseif r(2) == v, mxy = [ij(1)+1  r(2)*dy(2)+xy(3)];
        elseif r(3) == v, mxy = [xy(1)+r(3)*dy(1)    ij(2)];
        else              mxy = [xy(1)+r(4)*dy(1)  1+ij(2)]; 
        end
        
        img(ij(2),ij(1)) = img(ij(2),ij(1)) + sqrt(sum((xy(1,:)-mxy).^2)); 
        xy(1,:) = mxy;
    end
end

img = img ./ units(1); 



function show_density(rdat, anat, dd_img, xy_zero)

anat.node(:,1:2) = anat.node(:,1:2) + xy_zero;

cla
tools.loadAnatomy(anat, '-plot');
view([0 0 1]), tidyPlotForIllustrator % convert from 3D 
h = get(gca,'children'); 
set(h,'Color',[.1 .1 .1 .7],'LineWidth',1.5,'MArkerFaceColor',[.1 .1 .1])
set(gca,'children',flipud(get(gca,'children')))

hold on, axis image
imagesc(rdat.range,rdat.range,dd_img)
set(gca,'children',flipud(get(gca,'children')))

% Calibration image - check rounding of pixel values 

% xy_zero = [-560 0]
% dend_xy = bsxfun(@plus,anat.dendrites,xy_zero)*units(1) + units(2) + 0.5;
% p_x = (1:length(img)) + 0.5;

% clf
% subplot(1,2,1)
% imagesc(p_x,p_x,rdat.image), hold on, axis image ij
% plot(dend_xy(:,1),dend_xy(:,2),'.r-')

% subplot(1,2,2)
% imagesc(rdat.range,rdat.range,rdat.image), hold on, axis image ij
% plot(anat.dendrites(:,1)+xy_zero(1),anat.dendrites(:,2)+xy_zero(2),'.r-')

return




%% Anatomy Mutual information Analysis (imported from old masterfile)
function D = align_anatomy_RF(D) % HEKA  analysis (old FIG_3) 

do_manual = any(strncmpi(D,'-man',4)); % Manually do first pass of alignment
do_PDF = ~any(strncmpi(D,'-pdf',4)); % Save a PDF of the processing (default)
set_nK = any(strncmpi(D,'-nK',3)); % Default: 3-component
do_pmkmp = any(strncmpi(D,'-pmkmp',3)); % Don't use individual component colors

if do_PDF
  plot_folder = ['radon' filesep];
  if exist([tempdir plot_folder], 'file'), rmdir([tempdir plot_folder], 's'), end
  mkdir([tempdir plot_folder])
  print_options = {'-move','-clps'}; % '-portrait', '-resize'
end

if set_nK
    nK = D{strncmpi(D,'-nK',3)}; 
    nK = str2double(regexp(nK,'\d+','match','once'));
    if isnan(nK), nK = 2; end
else nK = 3;
end

physiol_list = HEKA_file_list(D); 

T = table; 

if false
    
    T.file = physiol_list(10); 
else
    T.file = physiol_list(:); 
end


%%
new_column = cell(size(T.file)); 

T.XY = new_column;
T.best_MI_id = zeros(size(new_column)); 
T.MI_img = new_column;
T.Hm_img = new_column;
T.MI_img_XY = new_column;
T.iterative_MI = new_column;
T.iterative_Hm = new_column;
T.fbp_MI = new_column;
T.fbp_Hm = new_column;
T.blur_radius = zeros(size(new_column));
T.radon_barwidth = zeros(size(new_column));
T.radon_tf = zeros(size(new_column));

T.Properties.UserData.Generated_At = datestr(now()); 

units.XY = {'X, µm','Y, µm'};
units.best_MI_id = 'component ID [psth.nnmf(nK=3)]'; 
units.MI_img = 'Mutual Information (as %) vs displacement';
units.Hm_img = 'image entropy vs displacement';
units.MI_img_XY = 'displacement, µm (X,Y)';
units.iter_MI = {'Mutual Information (as %)', 'at XY', 'max', 'imreg max', '... for each component'};
units.iter_Hm = {'image entropy', 'at XY', 'max', 'imreg max', '... for each component'};
units.fbp_MI = 'Mutual Information (as %) for FBP for each component';
units.fbp_Hm = 'image entropy (as %) for FBP for each component';

clear anat_list physiol_list new_column ii

%%

for ff = 1:height(T)

    close all
    if ~exist('ff','var'), ff = 15; end %#ok<FXSET> - '8.6.1 #13'
    Tools.loadPhysiology(T.file{ff},'-nnmf','-psth','-nK',nK,'-anatomy')
    if ~isfield(anatomy,'dendrites'), continue, end
    %%
    Tools.standardFigure('Landscape','Name','Anatomy MI')
    old_rdat = Tools.plot_radon_IMG(coeff,'ALG','OLD_radon');
    rdat = Tools.plot_radon_IMG(coeff,score,'ALG','sart'); %#ok<NODEF>
    if all(isnan(rdat.y_base)), continue, end, pause(0.01)
    
    %% Step 1: Get in the ballpark using imregconfig    
    figure(2), clf
    Tools.makeDendriticDensityImg([]) % Reset    
    img = Tools.makeDendriticDensityImg(rdat,anatomy,[0 0],'-plot');    
    
    h = get(gca,'Children');
    roi = [min(h(end-1).XData) max(h(end-1).XData) ...
           min(h(end-1).YData) max(h(end-1).YData)];
    roi = roi - [1 -1 1 -1]*range(xlim)/50; 
    axis(roi), axis off, pause(0.01)
        
    X = rdat.range; 
    
    best.id = []; 
    best.val = 0; 
    best.xy = [];
    
    iter_MI = zeros(nK,4); 
    iter_Hm = zeros(nK,4); 
    printInfo;
    
    for kk = 1:nK % For each component 
      %%
      printInfo('Computing image registration, k%d...\n', kk)
      rdat.image = rdat.images{kk};

      if do_manual
      
          figure(1), pause(0.01)
          h = flipud(get(gcf,'Children'));
          h = h(kk*3); 
          set(h,'Visible','on','XColor','r','YColor','r','LineWidth',2,'XTick',[],'YTick',[])
          hold(h,'on'), 
          if exist('xy_t','var')
               plot(xy_t(1),xy_t(2),'r+','Parent',h)
          else plot(0,0,'r+','Parent',h)
          end
          [xy_t(1),xy_t(2)] = ginput(1); 
          set(h.Children(1),'marker','o','xdata',xy_t(1),'ydata',xy_t(2))          
          axis(h,'off'), pause(0.01)
      else
          [optimizer, metric] = imregconfig('multimodal'); 

          for refinement = 1:5
            lastwarn('Reset Warnings to pick up images:regmex:registrationOutBoundsTermination', ...
                   'heka:okToContinue')
            xform = imregtform(img, rdat.image, 'translation', optimizer, metric);
            optimizer.InitialRadius = optimizer.InitialRadius / 2;

            [~,warning_message] = lastwarn; 
            if strcmp(warning_message,'heka:okToContinue'), break, end
          end
          xy_t = xform.T(3,[1 2]) * mean(diff(X)); 
      end
      
      [mi0,ih0] = Tools.mutualInformation(rdat.image, ...
                  Tools.makeDendriticDensityImg(rdat,anatomy,xy_t)); 
              
      get_MI_dxy = @(xy) -Tools.mutualInformation(rdat.image, ...
                           Tools.makeDendriticDensityImg(rdat,anatomy,xy));
      d_xy = fminsearch(get_MI_dxy, xy_t); 
      
      % Estimate point-spread function 
      xy_gauss = meshgrid( (-10:10) * mean(diff(X)));
      xy_gauss = sqrt(xy_gauss.^2 + xy_gauss' .^2);
      gauss = @(r) exp(-(xy_gauss/r).^2);
      gauss_smooth = @(r,im) conv2(im,gauss(r),'same');
      
      im_fixed = Tools.makeDendriticDensityImg(rdat,anatomy,d_xy);      
      get_MI_rb = @(r) -Tools.mutualInformation(rdat.image, ...
                                   gauss_smooth(r,im_fixed));
      r_b = fminbnd(get_MI_rb, 0, 40); 
      
      % Based on https://doi.org/10.3389/fnana.2017.00092, mouse bipolar
      % cells have axon arbors between 5 and 19 µm, so I'm happy to go up
      % to double that (no more)
      
      if 0 % show Gaussian smoothing examples
        %%
        clf, subplot(2,2,1) %#ok<UNRCH>
        im_fixed = Tools.makeDendriticDensityImg(rdat,anatomy,d_xy,'-plot');
        h = get(gca,'Children');
        roi = [min(h(end-1).XData) max(h(end-1).XData) ...
               min(h(end-1).YData) max(h(end-1).YData)];
        roi = roi - [1 -1 1 -1]*range(xlim)/50; 
        axis(roi), axis off, pause(0.01)
        
        g_radii = [0.3 1 3] * r_b; 
        for rr = 1:3
          subplot(3,2,2*rr), imagesc(X,X,gauss_smooth(g_radii(rr),im_fixed))
          axis image off, axis(roi)
        end
        
        subplot(2,2,3), imagesc(X,X,rdat.image), axis image off, axis(roi)
      end
            
      [mi,ih] = Tools.mutualInformation(rdat.image, gauss_smooth(r_b,im_fixed));       
      
      if (kk == 1) % mi/ih > best.val, 
        best.id = kk; 
        best.val = mi/ih;
        best.xy = d_xy; 
        best.rad = r_b; 
      end
       
      iter_MI(kk,2) = mi/ih;
      iter_Hm(kk,2) = ih;
       
      iter_MI(kk,3) = mi0/ih0;
      iter_Hm(kk,3) = ih0;
    end
    
    clear xy_t xform optimizer metric mi0 ih0 get_MI_d_xy roi 
    
    kk = best.id;
    rdat.image = rdat.images{kk};
    d_xy = best.xy; 
    
    printInfo('Best image registration: k%d\n', kk)
    
    %% Having picked a component, ...
    
    gx = (-12:12) * 5; 
   [gx,gy] = meshgrid(gx); 
    gx = gx + d_xy(1); 
    gy = gy + d_xy(2); 

    MI_surface = zeros(size(gx));
    Hm_surface = zeros(size(gx));

    Tools.mutualInformation(rdat.image, ... % reset internals
    Tools.makeDendriticDensityImg(rdat,anatomy,d_xy)); 
    
analysis.dendriticDensity

    for ii = 1:numel(gx) % compute mutual information 

      % break % Skip this calculation
      
      printInfo('Calculating MI surface [%0.1f%%]', 100*ii/numel(gx)) %%#ok<UNRCH>
      [mi,ih] = Tools.mutualInformation(rdat.image, gauss_smooth( best.rad,  ...
                Tools.makeDendriticDensityImg(rdat,anatomy,[gx(ii) gy(ii)])), true); 
      MI_surface(ii) = mi/ih;
      Hm_surface(ii) = ih; 
    end   

    if 0
        %% Plot computed surface 
        clf, hold on %#ok<UNRCH>
        imagesc(gx(1,:),gy(:,1),MI_surface * 100)
        colorbar, caxis([0 max(caxis)]), colormap pmkmp, axis off
        [val,ii] = max(MI_surface(:)); 
        
        plot(gx(ii),gy(ii),'.','Color',[.3 .3 .3],'MArkerSize',12)
        text(gx(ii)+3,gy(ii),sprintf('%0.2f%%',100*val),'Color',[.3 .3 .3])

        plot([0 50],max(ylim)+[2 2],'-','Color',[.3 .3 .3],'LineWidth',2,'Clipping','off') 
        text(25,max(ylim)+3,'50 µm', 'Color',[.3 .3 .3],'HorizontalAlignment','center','FontSize',12)
    end
    
    %%
    
    for kk = 1:nK % Get MI values for each component at best XY
        
      img = Tools.makeDendriticDensityImg(rdat,anatomy,d_xy); 
      [mi,ih] = Tools.mutualInformation(rdat.images{kk}, img);
                    
      iter_MI(kk,1) = mi/ih;
      iter_Hm(kk,1) = ih;
      
      [mi,ih] = Tools.mutualInformation(old_rdat.images{kk}, img);
                          
      iter_MI(kk,4) = mi/ih;
      iter_Hm(kk,4) = ih;
    end
    
    %% ... and save data to table 
    if exist('T','var')
        T.XY{ff} = d_xy;
        T.best_MI_id(ff) = best.id;
        T.MI_img{ff} = MI_surface;
        T.Hm_img{ff} = Hm_surface;
        T.MI_img_XY{ff} = [gx(1,:)' gy(:,1)];
        T.iterative_MI{ff} = iter_MI(:,1:3);
        T.iterative_Hm{ff} = iter_Hm(:,1:3);
        T.fbp_MI{ff} = iter_MI(:,4);
        T.fbp_Hm{ff} = iter_Hm(:,4);
        T.blur_radius(ff) = best.rad; 
        
        % Also get the key stimulus parameters:         
        expo_param = @(a,b) cellfun(@(e)e.Data{a}(b),expoData.passes.events(2:2:end));

        barY_pix = expo_param(3,5); 
        barY_deg = expo_param(4,5);
        deg_to_pix = median(barY_pix(barY_deg ~= 0) ./ barY_deg(barY_deg ~= 0)); 
        um_per_pixel = (2.342818332); % um per pixel

        T.radon_barwidth(ff) = median(expo_param(4,3) * deg_to_pix * um_per_pixel); 
        T.radon_tf(ff) = round(median(expo_param(1,3))*10)/10;

        save('HEKA MI analysis (anatomy).mat','T','units')
    end
    
    %%
    % de_units = (1:length(rdat.range))  / [ rdat.range; 0*rdat.range+1 ];    
    
%     dend_xy = bsxfun(@plus,anatomy.dendrites,d_xy);    
%     [px,py] = meshgrid(rdat.range);
%     
%     g_dist = arrayfun(@(x,y) nanmin(sum(bsxfun(@minus,[x y],dend_xy).^2,2)), px, py);
%     
%     
%     figure, semilogx(g_dist(rdat.mask),rdat.images{2}(rdat.mask),'.')

    %% Make PDF page (PNAS figure 3)
    txt_style = {'Color',[.3 .3 .3],'HorizontalAlignment','center','FontSize',9,'VerticalAlignment','top'};

    if exist('MI_surface','var')
        [~,ii] = max(MI_surface(:)); 
        d_xy = [gx(ii) gy(ii)];        
    end
    img = Tools.makeDendriticDensityImg(rdat,anatomy,d_xy);  
    dA = mean(diff(rdat.range)).^2;  % µm² / pixel
        
    img = gauss_smooth(best.rad, img / dA) / sum(sum(gauss(best.rad))); % now µm / µm²
    map = rdat.image / dA * max(score(:,1)); %%#ok<NODEF> % imp/s/µm²
    the = @(n,r) anatomy.(n)(:,r) + d_xy(r);

    close all
    Tools.standardFigure('Name','Anatomy RADON','Landscape','Tools')
    roi = [gx(ii)+anatomy.X([1 end]) gy(ii)-anatomy.Y([1 end])];

    subplot(3,4,1)
    image(anatomy.X, anatomy.Y, (anatomy.img + 0.25) / 1.25)
    axis image xy off, hold on, axis(axis)
    plot([0 -50],min(ylim)-[6 6],'-','Color',G(3),'LineWidth',2,'Clipping','off') 
    text(-25,min(ylim)-6,'50 µm', txt_style{:})

    % Traced anatomy
    subplot(3,4,5), hold on
    plot(the('dendrites',1), the('dendrites',2), 'Color',G(2),'LineWidth',0.8)
    if isfield(anatomy,'soma')
       fill(the('soma',1), the('soma',2),G(4),'EdgeColor',G(2))
    end
    plot(roi([1 2 2 1 1]),roi([3 3 4 4 3]),'-k')
    axis image ij off, hold on, axis(roi)   
    plot([0 -50],max(ylim)+[6 6],'-','Color',G(3),'LineWidth',2,'Clipping','off') 
    text(-25,max(ylim)+6,'50 µm', txt_style{:})
    
    subplot(3,4,[3 8])
    imagesc(rdat.range, rdat.range, map)
    axis image off, hold on, clim = caxis; 

    if do_pmkmp
    
        % "Component" color-map 
        c = pmkmp(201,'edge');
        caxis([-1 1]*max(abs(caxis)))
    else
        % "Component" color-map 
        c = interp1([-1;0.5;1], [W(1,5); C(1,:); C(1,:)/2], linspace(-1,1,64));
        c(isnan(c)) = 0.4; % "extrap" hit the boundry
    end
    colormap(gca,c);  axis(axis)
    plot(roi([1 2 2 1 1]),roi([3 3 4 4 3]),'-k')
    plot([0 100],max(ylim)+[10 10],'-','Color',[.3 .3 .3],'LineWidth',2,'Clipping','off') 
    text(50,max(ylim)+12,'100 µm', txt_style{:})
    
    
    
    
    % PSTH axes
    subplot(3,4,[11 12])
    bar(psth.time, score(:,best.id),1,'FaceColor',G(7),'EdgeColor','none'), axis(axis), hold on
    for ss = 1:nStimuli,
        rectangle('Position',stim_bar(ss,-0.1),'FaceColor',[0 0 0 0.5], 'EdgeColor','none')
    end
    tidyPlotForIllustrator, axis tight
    ylabel('PSTH component, imp/s')
    set(gca,'TickLength',[1 1]/80)

    % Plot right half ~ inserts and placeholder scatter graph
    sx = rdat.range >= roi(1) & rdat.range <= roi(2);
    sx = conv(double(sx),[1 1 1],'same')>0; 
    sy = rdat.range >= roi(3) & rdat.range <= roi(4);
    sy = conv(double(sy),[1 1 1],'same')>0;     
   
    % Plot inset of RF map
    subplot(3,4,2),
    imagesc(rdat.range(sx),rdat.range(sy), map(sy,sx))
    caxis(clim), axis image off, axis(roi), hold on, colormap(c)
    
    p = get(gca,'Position'); 
    ch = colorbar; ch.Position = p + [16.3 2.3 -15 -4.5]/100;
    plot([0 -50],max(ylim)+[6 6],'-','Color',[.3 .3 .3],'LineWidth',2,'Clipping','off') 
    text(-25,max(ylim)+6,'50 µm', txt_style{:})
    
    if do_pmkmp
       cl = caxis; 
       ch.Ticks = unique([cl ch.Ticks]);
       caxis([-1 1]*max(abs(cl)))        
    end
    
    % Plot dendritic density
    subplot(3,4,6)
    imagesc(rdat.range(sx),rdat.range(sy), img(sy,sx))
    axis image off, axis(roi), hold on, % colormap(gca,bone)
    plot(the('dendrites',1), the('dendrites',2), 'Color',[0 0 0 .3],'LineWidth',0.5)
    if isfield(anatomy','soma')
        fill(the('soma',1), the('soma',2),  [.4 .4 .4], 'EdgeColor',[.2 .2 .2])
    end
    plot([0 -50],max(ylim)+[6 6],'-','Color',[.3 .3 .3],'LineWidth',2,'Clipping','off') 
    text(-25,max(ylim)+6,'50 µm', txt_style{:})
    
    px = linspace(0,2*pi,65); 
    py = sin(px)*best.rad; px = cos(px)*best.rad;    
    plot(px + roi(1),py+roi(3),'-','Color',G(3),'LineWidth',1.2,'Clipping','off')
    
    c = 1-(summer(64).^2); 
    w = linspace(1,0,64)'.^2  * [3 3 3];
    colormap(gca, (c + w)./(1+w))
    
    p = get(gca,'Position'); 
    ch = colorbar; ch.Position = p + [16.3 2.3 -15 -4.5]/100;
    plot([0 -50],max(ylim)+[6 6],'-','Color',[.3 .3 .3],'LineWidth',2,'Clipping','off') 
    text(-25,max(ylim)+6,'50 µm', txt_style{:})
    
    % % Plot scatterplot Y-data inset: naive FBP
    % subplot(2,4,3), tweak(gca,[0 -1 2 0]); 
    % imagesc(rdat.range(sx),rdat.range(sy), old_rdat.image(sy,sx))
    % caxis([nanmin(old_rdat.image(:)) nanmax(old_rdat.image(:))])
    % axis image off, axis(roi), hold on, 
    % plot(the('dendrites',1), the('dendrites',2), 'Color',[0 0 0 .3],'LineWidth',0.8)
    % fill(the('soma',1), the('soma',2),  [.4 .4 .4], 'EdgeColor',[.2 .2 .2])

    % Scatterplot axes
    
    kz = quantile(img(img>0),0.05);
    subplot(3,4,[9 10]), hold on, % tweak(gca,[-5 -18 8 16])
    plot(-kz*(1+rand(sum(img(:) == 0),1)), map(img == 0), '.','Color',G(8))
    plot(img(img>0), map(img > 0), '.','Color',G(5))
    axis tight, tidyPlotForIllustrator, hold on    
        
    [r,p] = corr(img(img>0),map(img>0));
    p = min(1,p*height(T)); % Bonferroni-correction
    px = [min(xlim)/2 max(xlim)];
    r_line = [img(img>0) img(img>0)*0+1] \ map(img > 0);
    plot(px, [px; 1 1]*r_line,'-','Color',[0 0 0 0.5],'LineWidth',1.1)
    text(min(xlim)/2,max(ylim),sprintf('r = %0.3f (p = %0.4f)',r,p),'VerticalAlignment','top')
    
    xlabel('Dendrite density, \mum/\mum^2')
    ylabel('Receptive field, imp/s/\mum^2')

%     % MI surface axes
%     subplot(6,3,10), cla, hold on
%     imagesc(-gx(1,:),gy(:,1),MI_surface * 100)
%     colormap(gca,summer.^1.2), axis tight image off ij
%     plot(-gx(ii),gy(ii),'.','Color',[.3 .3 .3],'MArkerSize',12)
%     text(8-gx(ii),gy(ii),sprintf('%0.2f%%',100*val),'Color',[.3 .3 .3])
% 
%     plot([0 50],max(ylim)+[6 6],'-','Color',[.3 .3 .3],'LineWidth',2,'Clipping','off') 
%     text(25,max(ylim)+12,'50 µm', 'Color',[.3 .3 .3],'HorizontalAlignment','center','FontSize',9)

    suptitle(strrep(T.file{ff},'_','\_'))
        
    %% ... and save page to PDF
    if ~do_PDF, continue, end
    set(gcf,'PaperPositionMode','auto');
    output_file = sprintf('%s%03d_f03_plot.ps', plot_folder, ff);
    figs_to_ps_2015(gcf,output_file,print_options{:});
    pause(1)
    
end

% Remove skipped entries (no / missing anat)
T(T.best_MI_id == 0,:) = [];
save('HEKA MI analysis (anatomy).mat','T','units')

%%

% cd('..\LGN Radon\')
% now for the big finale, combine all ps files into a single pdf
if do_PDF, 
    tic, clc
    combine_PStoPDF([tempdir plot_folder], 'HEKA Mutual Information analysis.pdf')
    clc, toc, 
end

