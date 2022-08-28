function dat = plot_radon_IMG(dat,varargin)

named = @(n) strncmpi(varargin,n,length(n));
do_profiles = false; 

if ~isstruct(dat)

    if nargin > 1 && isnumeric(varargin{1})
        
        do_profiles = true; 
        profile = varargin{1};
        
        time = evalin('caller','time');
        if length(time) ~= size(profile,1)
            time = evalin('caller','psth.time'); 
            
            % convert to PSTH bins
            time = [reshape(time,1,[]) - mean(diff(time))/2 ...
                    time(end) + mean(diff(time))/2];
            time = reshape([1;1] * time,1,[]);
            profile = cat(1, permute(profile,[3 1 2]), permute(profile,[3 1 2]));
            profile = reshape(profile,[],size(profile,3));
            profile = [0*profile(1,:); profile; 0*profile(1,:)];
        end
    end
    
    wave = dat;
    expoData = evalin('caller','expoData');
    dat = struct; 
    
    % Extract Stimulus information
    A = cellfun(@(p) p.Data{3}(1), expoData.passes.events(2:2:end));
    X = cellfun(@(p) p.Data{3}(2), expoData.passes.events(2:2:end));
    X(A >= 180) = -X(A >= 180); 
    A(A >= 180) =  A(A >= 180) - 180; 

    [dat.ori,~,idx] = unique([A' X'],'rows');
    
    if any(named('-lgn')), um_per_pixel = 1; % degrees of visual field
    else                   um_per_pixel = 2.342818332; % Screen conversion
    end

    dat.x   = dat.ori(:,2) * um_per_pixel; % um per pixel
    dat.ori = dat.ori(:,1);

    for ii = 1:length(dat.ori)    
        dat.y_all(ii,:) = mean(wave(idx == ii,:), 1);
    end
    
    clear expoData wave
end





%%
fprintf('Plotting Radon ...\n')

pause(0.1), clf, nK = size(dat.y_all,2);

if do_profiles
  subplot(1,4,1); p1 = get(gca,'Position')./[1 1 1 nK]; delete(gca);
  subplot(1,4,[2 3]); p2 = get(gca,'Position')./[1 1 1 nK]; delete(gca);
  subplot(1,4,4); p3 = get(gca,'Position')./[1 1 1 nK]; delete(gca);
else
  subplot(1,2,1); p2 = get(gca,'Position')./[1 1 1 nK]; delete(gca);
  subplot(1,2,2); p3 = get(gca,'Position')./[1 1 1 nK]; delete(gca);
end

% cl = lines(nK);
nX = length(unique(dat.x));
    
dat.images = cell(size(dat.y_all(1,:)));
dat.y_base = nanmedian(dat.y_all,1);

for kk = 1:nK    
    %% Axis 1 - Profile
    if do_profiles
      axes('Position',p1+[0 (nK - kk)*p1(4) 0 0])    
      plot(time,profile(:,kk),'Color',[.1 .1 .1],'LineWidth',1.2)
      ylim(ylim);
     
%       for ss = 1:dat.stim.count,
%         rectangle('Position',dat.stim.bar(ss-1),'FaceColor',[.5 .5 .5 .4], 'EdgeColor','none')
%         rectangle('Position',dat.stim.bar(ss-0.5),'FaceColor',[0  0  0  .4], 'EdgeColor','none')
%       end   

      xlim([min(time) max(time)] * [1.01 -0.01; -0.01 1.01])

      if kk == nK, set(gca,'XTickLabel', strcat(get(gca,'XTickLabel'),' s'))
      else         set(gca,'XTickLabel',{})
      end

      yl = ylim; h = gca;
      if time(1) == time(2), y_unit = ' imp'; else y_unit = ' mV'; end

      h.YTick = unique([yl(1) 0 yl(2)]);
      text(xlim*[1.02;-.02],ylim*[.98;.02],[h.YTickLabel{1} y_unit],'VerticalAlignment','Bottom','HorizontalAlignment','right')
      text(xlim*[1.02;-.02],ylim*[.02;.98],[h.YTickLabel{end} y_unit],'VerticalAlignment','Top','HorizontalAlignment','right')    
      h.YTickLabel = []; grid on
      set(gca,'UserData',kk)
    end

    %% Axis 2 - Coefficients
    axes('Position',p2+[0 (nK - kk)*p2(4) 0 0])
    imagesc(dat.x(1:nX), dat.ori(1:nX:end),reshape(dat.y_all(:,kk),nX,[])')
    colorbar
    
    if kk == nK, set(gca,'XTickLabel', strcat(get(gca,'XTickLabel'),'°'))
    else         set(gca,'XTickLabel',{})
    end
    
    set(gca,'YTick',dat.ori(1:nX:end),'TickLength',[0 0])
    set(gca,'YTickLabel',strcat(get(gca,'YTickLabel'),' µm'))
    set(gca,'UserData',kk)
    
    %% Axis 3 - Radon Transform
    axes('Position',p3+[0 (nK - kk)*p2(4) 0 0])
    
    if any(named('ALG'))        
        dat.algorithm = varargin{find(named('ALG'))+1};
    end
    
    dat.y = dat.y_all(:,kk) - dat.y_base(kk);
    dat = Tools.iRadon(dat);

    imagesc(dat.range,dat.range,dat.image)
    axis image off xy
    set(gca,'UserData',kk)

    dat.images{kk} = dat.image;
    if isfield(dat,'image_0')
        dat.img_FBP{kk} = dat.image_0;
    end
end

x = max(abs(dat.x));

hold on, axis(axis)
plot([0 100],-1.05*[x x],'-','LineWidth',1.5,'Color',[.3 .3 .3],'Clipping','off')
text(50,-1.15*x,sprintf('100 µm'),'FontSize',13,'HorizontalAlignment','center','Color',[.3 .3 .3])

% %% Latent axes - explained 
% axes('Position',p3+[0 (nK - kk-1)*p2(4)+0.02 0 -0.035])
% latent = dat.latent; explained = sum(latent(1:nK))/sum(latent)*100;
% 
% plot(1:numel(latent), cumsum(latent)/sum(latent)*100, 'LineWidth',1.3)
% hold on, plot([nK nK nan xlim],[ylim nan explained*[1 1]],'LineWidth',1.2)
% xlim([1 numel(latent)]), title('\rmVariance Explained')
% set(gca,'YTickMode','manual','YTickLabel',strcat(get(gca,'YTickLabel'),'%'),'FontSize',7)

if nargout == 0, clear, end

end
