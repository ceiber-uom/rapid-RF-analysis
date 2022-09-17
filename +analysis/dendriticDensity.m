function img = dendriticDensity(rdat,anat,xy_zero,do_plot)

if nargin < 3, xy_zero = [0 0]; end
if nargin < 2, anat = evalin('caller','anatomy'); end
if nargin < 1, rdat = evalin('caller','rdat'); end

% this is going to be called in a loop ... precompute and optimise
persistent units dd_indices 
if isempty(rdat), units = []; return, end % reset
if isempty(units)
    units = (1:length(rdat.range))  / [ rdat.range; 0*rdat.range+1 ];    

    dd_indices = find(~( isnan(anat.dendrites(1:end-1,1)) | ...
                         isnan(anat.dendrites(2:end,1)))); 
    dd_indices = reshape(dd_indices,1,[]);
end

if ~isfield(anat,'dendrites'), 
    error('in order to run "%s", traced dendritic anatomy is required.', mfilename)
end

%%
dend_xy = bsxfun(@plus,anat.dendrites,xy_zero)*units(1) + units(2) + 0.5;
% pixelized co-ordinates: a point is in px(x=1) if 1 < p.x < 2

img = zeros(size(rdat.image)); 
for dd = dd_indices
    
    xy = dend_xy(dd+[0 1],:);
    
    while true
        
        dy = diff(xy);
        ij = floor(xy(1,:));        
        
        if any(ij < 1 | ij([2 1]) > size(img)), break, end
        
        r = (ij([1 1 2 2])-xy([1 1 3 3])+[0 1 0 1])./dy([1 1 2 2]); 
        v = min(r(r>0)); 
        
        if v > 1, 
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

%%
if nargout > 0 && (nargin < 4 || ~any(do_plot)), return, end

%%

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



%%

cla
the = @(n,r) anat.(n)(:,r) + xy_zero(r);
imagesc(rdat.range,rdat.range,img), hold on, axis image ij
plot(the('dendrites',1), the('dendrites',2), 'Color',[0 .5 0 .8],'LineWidth',1)
if isfield(anat,'soma')
  fill(the('soma',1), the('soma',2),  [.6 .8 .6], 'EdgeColor',[0 .5 0])
end

c = 1-(summer.^3); 
w = linspace(1,0,64)'.^3  * [4 4 4];
colormap(gca, (c + w)./(1+w))
