function plot_anatomy(date, rdat, anat, xy, gm, varargin)

f3 = figure;
ax_h(1) = subplot(1,2,1);
d = rdat.images{1};
% See notes: PRJ-VisionLab\Elissa_Belluccini\Notes\Orienting Cell Morphology and RF Map.docx
% ASB correction: Rotate image 90 deg anticlockwise and flip about horizontal AND vertical axis
if (date < 20210526) 
    map = imagesc(rdat.range,rdat.range,fliplr(flipud(imrotate(d,-90))));
else
% MFB correction: Rotate image 90 deg anticlockwise and flip about horizontal axis
    map = imagesc(rdat.range,rdat.range,flipud(imrotate(d,-90)));       
end

% map = imagesc(rdat.range,rdat.range,d);
axis square
axis image off xy
hold on
if ~isempty(varargin)       
    midx = varargin(1);
    midy = varargin(2);
else
    disp('Select receptive field centre');
    [midx,midy] = ginput(1);
end
% set(map.Parent,'CLim',[-1 1] * max(abs([map.Parent.CLim])))
% cm = interp1((-5:5)', redbluecmap, linspace(-5,5,101)); 
% ca = gca;
% colormap(ca,cm)
c = colorbar; 
c.TickDirection = 'out';
c.Box = 'off'; 
c.Location = 'westoutside';
cpos = c.Position;
c.Position = cpos.*[0.3,0.75,1.2,1.82]; % hard-coded 
axis image on xy
tidyPlotForIllustrator

ac_sz = size(anat.CData);
xdif = mean(arrayfun(@(xx) map.XData(xx)-map.XData(xx-1), 2:size(map.XData,2)));
ydif = mean(arrayfun(@(yy) map.YData(yy)-map.YData(yy-1), 2:size(map.YData,2)));
map_um_per_pix = mean([xdif,ydif]);

anat_width_um = abs(xy(1))+abs(xy(2));
anat_height_um = abs(xy(3))+abs(xy(4));
xdif = mean(arrayfun(@(xx) anat.XData(xx)-anat.XData(xx-1), 2:size(anat.XData,2)));
ydif = mean(arrayfun(@(yy) anat.YData(yy)-anat.YData(yy-1), 2:size(anat.YData,2)));
anat_um_per_pix = mean([abs(xdif),abs(ydif)]);

roi_map_w_pix = anat_width_um/map_um_per_pix;
roi_map_h_pix = anat_height_um/map_um_per_pix;

% Zoom in on RF
% xy2 = xy + [-10,10,10,-10];

% For 20200116_Cell_02
% xy2 = xy2 + [0,0,-50,50];

[min1,i1] = min(abs(rdat.range-xy2(1))); %i1 = i1-1;
[min2,i2] = min(abs(rdat.range-xy2(2))); %i2 = i2+1;
[min3,i3] = min(abs(rdat.range-xy2(3))); %i3 = i3+1;
[min4,i4] = min(abs(rdat.range-xy2(4))); %i4 = i4-1;

ax_h(2) = subplot(1,2,2);
m = imagesc(rdat.range(i1:i2),rdat.range(i4:i3),map.CData(i4:i3,i1:i2));
axis image off xy
axis image on xy 
hold on

imagesc(linspace(xy(1),xy(2),ac_sz(2)),linspace(xy(3),xy(4),ac_sz(1)),anat.CData,...
    'AlphaData',anat.AlphaData);

C = lines(7);
oneSD_circle = [cos(linspace(0,2*pi,61));sin(linspace(0,2*pi,61))]'; 
style = {'Color',C(1,:),'LineWidth',1.2};
xy = gm.center_xy(1,:) + gm.gauss_radius(1) * oneSD_circle;
xy = [midx{1},midy{1}] + gm.gauss_radius(1) * oneSD_circle;
plot(xy(:,1),xy(:,2),'-',style{:})


end