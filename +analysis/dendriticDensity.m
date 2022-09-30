
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
