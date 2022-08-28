function [img,xy] = loadAnatomy(expoData)
% This also functions as a utility to get the anatomy processed into a
% format which is suitable for other uses by MATLAB, if that hasn't been
% done already. Originally called "plotAnatomy" because plotting the
% anatomy is the default behaviour. 
% 
% Usage: [img,xy] = loadAnatomy(file_name) using "YYYYMMDD_Cell_X" returns
%                               the anatomy image and bounding box for the
%                               requested cell. If it isn't found, the user
%                               will be walked through the process of
%                               creating the needed files. 
%        [img,xy] = loadAnatomy grabs "expoData" from the calling workspace
%                               and uses expoData.FileName. 
% calling loadAnatomy(...) or [img] = loadAnatomy(...) causes this to
%     generate a translucent image of the cell anatomy in the current axes,
%     which should be on top of an existing radon RF map, using the data in
%     "Anatomy.xlsx" to align the anatomy and physiology, and possibly
%     return a handle to the image object. Again, if the requested cell
%     isn't found, the user will be walked through the process of creating
%     the needed files. 
% calling loadAnatomy -plot will cause this to exit gracefully if the cell
%     is not found, preventing the program from hanging waiting for user
%     input. In this mode, "expoData.FileName" in the local workspace
%     specifies the cell identity. 


    mode_PLOT = false;
    if exist('expoData','var') && ischar(expoData)
        mode_PLOT = strncmpi(expoData,'-plot',5) || ...
                    strncmpi(expoData,'-check',3); % Backwards-compatible
        if ~mode_PLOT
            cmd = expoData;
            expoData = struct;
            expoData.FileName = cmd;
        else clear expoData
        end
    end

    if ~exist('expoData','var')
      if evalin('caller','exist(''expoData'',''var'')')
        expoData = evalin('caller','expoData');
      else 
        clear % Test mode
        load('.\MAT\20180703_Cell_1 #9[Radon_Flicker_ACH].mat')
        
        figure(1), clf
        Tools.plot_radon_IMG(hekaData.PassData(4000,:)') %%#ok<NODEF>
      end
    end
    
    %%    
    here = strrep(fileparts(mfilename('fullpath')),'+Tools','Anatomy');
    persistent cell_table screen_table
    if isempty(cell_table)
        cell_table = readtable([here filesep 'Anatomy.xlsx'],'Sheet',1);
        screen_table = readtable([here filesep 'Anatomy.xlsx'],'Sheet',2);
    end

    exp_date = str2double(regexp(expoData.FileName,'\d{8}','match','once')); 
    cell_no = str2double(regexp(expoData.FileName,'(?<=Cell_)\d+','match','once'));
    
    if isnan(cell_no), 
        cell_no = str2double(regexp(expoData.FileName,'\d+(?=#)','match','once'));
    end
    
% UNIX Time: seconds elapsed since 00:00:00 UTC, Thursday, 1 January 1970
%     time = hekaData.FileHeader.StartTime/24/60 %  + ...
%             datenum(1970,1,1); 

    %%
    jpg_um_per_px = (0.157728707); % Hard-coded
    sel = (screen_table.Date == exp_date);
    ref = find(sel & screen_table.Cell == cell_no,1); 
    sel = sel & screen_table.Cell == 0;
    
    if ~isempty(ref) % Interpolate between stimulus images
    
        tf = screen_table.Time(sel); 
        tr = screen_table.Time(ref);
        
        if isnan(tr), tr = nanmean(tf); end
        
        s1 = find(sel & screen_table.Time == max([tf(tf <= tr);0]),1); 
        s2 = find(sel & screen_table.Time == min([tf(tf >= tr);1]),1); 
        
        if isempty(s1), s1 = s2; end
        if isempty(s1), s1 = find(sel & (screen_table.Time == min(tf)),1); end
        if isempty(s2), s2 = s1; end
        
        tf = screen_table.Time([s1 s2]); 
        if s1 == s2, ij = 0; else ij = (tr-tf(1))/(tf(2)-tf(1)); end
        
        x0 = (screen_table.X_px(s1) *    ij   + ...
              screen_table.X_px(s2) * (1-ij)) * jpg_um_per_px; 

        y0 = (screen_table.Y_px(s1) *    ij   + ...
              screen_table.Y_px(s2) * (1-ij)) * jpg_um_per_px; 
          
          
        xr = (screen_table.X_px(ref)) * jpg_um_per_px; % Default: middle of screen
        yr = (screen_table.Y_px(ref)) * jpg_um_per_px; 

        clear tr tf s1 s2 ij        
    else
        sel = sel & screen_table.Cell == 0;
        
        if ~any(sel),
            x0 = []; y0 = []; % Kick this down till after we check to see if there's filled anatomy for this cell
        else
            sel = sel & (screen_table.Time == min(screen_table.Time(sel))); 
            x0 = screen_table.X_px(sel) * jpg_um_per_px;
            y0 = screen_table.Y_px(sel) * jpg_um_per_px;
        end
        
        xr = (1360 / 2) * jpg_um_per_px; % Default: middle of screen
        yr = (1024 / 2) * jpg_um_per_px; 
        
    end
    
    % Get match from table of cell anatomies
    sel = (cell_table.Date == exp_date & cell_table.Cell == cell_no);    
    
    if ~any(sel)
        warning('%d Cell %02d not found in Anatomy.xlsx', exp_date, cell_no)
        if mode_PLOT, clear, return, end
        [cell_table, sel] = manual_import(cell_table,exp_date,cell_no); 
    end

    if isempty(x0) ; % Anatomy but not screen coordinates found
        warning('%d not found in screen_table', exp_date),
        x0 = (1360 / 2) * jpg_um_per_px;
        y0 = (1024 / 2) * jpg_um_per_px; 
    end
    
    png_um_per_px = cell_table.Scale(sel); 
    xc = cell_table.X_Cell(sel) * png_um_per_px;
    yc = cell_table.Y_Cell(sel) * png_um_per_px;
    
    %%
    
    img = imread(fullfile(here,'Filled Anatomy',cell_table.File{sel}));
    
    xL = (xr - xc - x0); % Left edge of image in µm
    yT = (y0 + yc - yr); % Top edge of image in µm
    xR = xL + size(img,2) * png_um_per_px; % Right edge
    yB = yT - size(img,1) * png_um_per_px; % Bottom edge

    img = double(img) / quantile(double(img(:)), 0.98);     
    bg = quantile(img(img > 0),0.02); 
    img = (img-bg)/(1-bg);     
    img(img < 0) = 0;  
    img(img > 1) = 1;     
    xy = [xL xR yT yB]; 
    
    % pow = 1.5;
    
    if nargout == 2, return, end
    
    hold on, axis(axis)
    img = imagesc(linspace(xL,xR,size(img,2)), ...
            linspace(yT,yB,size(img,1)), img, ...
            'AlphaData', min(1,mean(img,3)).^1.5, ...
            'UserData', cell_table.File{sel});
    
    %%
    if nargout == 0, clear, end
    return
    %%  DEBUG display information

    disp('Debugging...') %#ok<UNRCH>
    
    plot([0 0 nan xlim],[ylim nan 0 0],'k-')    
    plot(xL,yT,'ro','MarkerFaceColor','r')
    plot([xL xR xR xL xL],[yT yT yB yB yT],'r-')

    plot(-x0,y0,'go','MarkerFaceColor','g')
    plot(-x0 + [0 0 1360 1360 0 nan 0 1360 nan 680 680] * jpg_um_per_px, ...
          y0 - [0 1024 1024 0 0 nan 512 512 nan 0 1024] * jpg_um_per_px, 'g-')
    
      
    %% Show reference images
    figure(2), clf, subplot(1,2,1)    
    [sz_fn,sz_fp] = uigetfile('*.jpg');    
    sz_img = imread([sz_fp sz_fn]);

    image(sz_img), axis image
    hold on, plot(x0 / jpg_um_per_px, y0 / jpg_um_per_px, 'ro'); 
    plot([680 680 nan 1 1360], [1 1024 nan 512 512], 'r-')
    
    subplot(1,2,2)
    [sz_fn,sz_fp] = uigetfile('*.png',[],sz_fp);    
    sz_img = imread([sz_fp sz_fn]);

    image(sz_img), axis image
    hold on, plot(x0 / jpg_um_per_px, y0 / jpg_um_per_px, 'ro'); 
    plot([680 680 nan 1 1360], [1 1024 nan 512 512], 'r-')
    
    
end



function [T, sel] = manual_import(T,exp_date,cell_no)
    
% clearvars -except fp
% if ~exist('fp','var') || ~ischar(fp),  end

f_path = '.\Anatomy\Filled Anatomy\';

do_automatic = exist('exp_date','var');

if do_automatic
      f_name = uigetfile({'*.png';'*.jpg'},'Select Image',f_path);
else [f_name, f_path] = uigetfile({'*.png';'*.jpg'},'Select Image',f_path);
    
    T = readtable('Anatomy\Anatomy.xlsx','Sheet',1);
    exp_date = str2double(regexp(f_name,'\d{8}','match','once')); 
    cell_no = str2double(regexp(f_name,'(?<=Cell_)\d+','match','once'));
end

img = imread([f_path f_name]);
img = mean(img,3);  ok = img >= quantile(img(:),0.75); 
img = img ./ quantile(img(ok),0.99);
img(img > 1) = 1; % img(~ok) = NaN;

if size(img, 1) > 2000  
    pow = 1.5; 
else % Assuming 1024 x 1360
    pow = 1;
end

clf, imagesc(img .^ pow), hold on, colormap bone, axis image

if any(T.Date == exp_date & T.Cell == cell_no)
    warning('%d_Cell_%02d found at row %d in Anatomy.xlsx',exp_date,cell_no,...
                            find(T.Date == exp_date & T.Cell == cell_no,1))
end

%%

clf, imagesc(img .^ pow), hold on, colormap bone, axis image

key = 0; 
while key ~= 10 

    if key == 0, title('Select 50 µm scalebar'), end
    
    [x,y,key] = ginput(1); 
    if isempty(key)
        if exist('pix_to_um','var'), break
        else continue
        end
    elseif key == 27, error('Cancelled by user')
    end
    
    
    x = round(x);
    y = round(y);
    ref = img(y,x); 

    ok = arrayfun(@(p) all(img(y,min([p x]):max([p x])) == ref), 1:size(img,2));
    x = [find(ok,1) find(ok,1,'last')];
    pix_to_um = 50 / diff(x); % sqrt(diff(x).^2 + diff(y).^2);     
    
    clf, imagesc(img .^ pow), hold on, colormap bone, axis image
    plot(x,[y y],'r-')

    title(sprintf('pix\\_to\\_um = %0.9f\nEnter to continue', pix_to_um))
    
end

clear ref ok x y key 


key = 0; 
while key ~= 10 
    if key == 0, title('Select Cell body'), end
    [x,y,key] = ginput(1); 
    if isempty(key), 
        if exist('xc','var'), break
        else continue
        end
    elseif key == 27, error('Cancelled by user')
    end
    
    xc = x; yc = y; 
    plot(x,y,'ro')
    text(x,y,sprintf('%0.1f\n%0.1f',x,y),'Color','r','FontSize',7)
    title(sprintf('pix\\_to\\_um = %0.9f\nEnter to continue', pix_to_um))    
end


sel = height(T) + 1; 

T(sel,{'Date'}) = {exp_date}; % Needed to extend table
T.Cell(sel) = cell_no; 
T.File{sel} = f_name; 
T.Scale(sel) = pix_to_um;
T.X_Cell(sel) = xc;
T.Y_Cell(sel) = yc;

writetable(T, 'Anatomy\Anatomy.xlsx','Sheet',1);

if ~do_automatic, 
    clear x y xc yc cell_no exp_date f_name key pow img sel pix_to_um ans do_automatic
end

end % import



%% Use this to make PNG images monochrome and do other cleanup
function deColor  %#ok<DEFNU>
%%

do_preview = ~exist('do_preview','var') || do_preview; %#ok<NODEF>

if ~exist('p','var'), p = './Anatomy/Filled Anatomy'; end
[f,p] = uigetfile('*.png',[],p);
img = imread([p f]);

if ~do_preview
  if all(img(1,1,:) == 255) % white diagonal outside, because AI doesn't
                          % set the black background correctly
    mv = quantile(img(:),0.85);
    for ii = 1:size(img,2);     
        jj = find(img(:,ii,1) < mv,1) + 2; 
        if isempty(jj) || jj >= size(img,1),
            img(:,ii,:) = 0;
            continue
        end
        
        img(1:jj,ii,:) = 0;
        jj = find(img(:,ii,1) < mv,1,'last') - 2; 
        img(jj-2:end,ii,:) = 0;         
    end
    
  end % if (255)

  img(:,:,3) = img(:,:,1); 
  img(:,:,2) = img(:,:,1);
end

% img = img(2:3:end,2:3:end,:);

image(img), axis image, title(f,'Interpreter','none')


%%
imwrite(img,[p f]); 
return

%%

% ok = any(img(:,:,1) ~= 255,1);
% img = img(ok,:,:);

% ok = any(img(:,:,1) ~= 255,2);
% img = img(:,ok,:);

% img = max(double(img)-46,0) / (255-46);


%% Clean up specific crap

[x,y] = ginput(1); %#ok<UNRCH>
x = round(x); y = round(y);
img(y:end,x,:) = 0;

% img(:,x,:) = 0;

h = get(gca,'Children');
set(h,'CData',img)

%% De-spackle

img = imread([p f]);

rr = 5; 
spc = zeros(size(img(:,:,1))); 

for ii = 1:size(img,2)
    idx = (ii-rr:ii+rr);
    idx(idx < 1) = [];
    idx(idx > size(img,2)) = [];     
    rows = mean(img(:,idx,1) > 100,2);
    spc(:,ii) = conv(rows,ones(2*rr+1,1)/(2*rr+1),'same');
end

lvl = 0.125; 
img(:,:,2) = min(spc / lvl, 1) * 255; % (spc > 0.1) * 255;

% c2 = conv2(img(:,:,1),ones(13)/13.^2,'same'); 

img(:,:,3) = cast(double(img(:,:,1)) .* min(spc / lvl, 1),'like',img);

imagesc(img(:,:,3)), axis image, colormap bone

%% 

img(:,:,2) = img(:,:,3);
img(:,:,1) = img(:,:,3); 
imagesc(img), axis image


end


