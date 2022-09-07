function dat = plot_radon_IMG(dat,varargin)
% plots.plot_radon_IMG( data, ... )
% plots.plot_radon_IMG( activations, [waves, time], ... )
% 
% Generates a standardised plot of the measured sinogram and and the 
% corresponding spatial map (generated using analysis.inverseRadon) 
% 
% if [waves, time] is also supplied the component waveforms coresponding to
% each spatial map. for PCA or NNMF components, the waveform is important
% for interpreting the spatial map. 
% 
% Options:
% -no-check : skip utils.prepareRadon (if using a data structure which is 
%                  already formatted correctly, can be more efficient)
% -no-base  : do not subtract baselines (median sinogram value) before
%                  generating spatial map. 
% -alg [algorithm] : specify algorithm for spatial map generation 
%                     (see analysis.inverseRadon for details)
% -latent [latent vector] : plot also fraction of variance explained 
%                     (see doc pca). [latent vector] is unnecessary if 
%                     data.latent supplied. 
% 
% v0.2 - 2 September 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(n) varargin{find(named(n))+1};

if ~any(named('-no-check')), dat = utils.prepareRadon(dat, varargin{:}); end

do_profiles = isfield(dat,'wave'); 


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
dat.y_base = nanmedian(dat.y_all,1); %#ok<NANMEDIAN> 
if any(named('-no-base')), dat.y_base(:) = 0; end

is_imp_s = (dat.time(1) == dat.time(2));

if is_imp_s, y_unit = 'imp'; 
elseif any(named('-units-V')), y_unit = 'V';
else y_unit = 'mV';
  if do_profiles && max(abs(dat.wave(:))) < 0.5
    dat.wave = 1e3*dat.wave; 
    disp('Converting wave units to mV (assumed: V)')
  end
end


for kk = 1:nK    
    %% Axis 1 - Profile
    if do_profiles
      axes('Position',p1+[0 (nK - kk)*p1(4) 0 0])    
      plot(dat.time,dat.wave(:,kk),'Color',[.1 .1 .1],'LineWidth',1.2)
      ylim(ylim);
     
%       for ss = 1:dat.stim.count,
%         rectangle('Position',dat.stim.bar(ss-1),'FaceColor',[.5 .5 .5 .4], 'EdgeColor','none')
%         rectangle('Position',dat.stim.bar(ss-0.5),'FaceColor',[0  0  0  .4], 'EdgeColor','none')
%       end   

      xl = [min(dat.time) max(dat.time)] * [1.01 -0.01; -0.01 1.01];       
      xlim(xl) % ±1% from minimum & maximum time 

      if kk == nK, set(gca,'XTickLabel', strcat(get(gca,'XTickLabel'),' s'))
      else         set(gca,'XTickLabel',{})
      end

      yl = ylim; h = gca;

      h.YTick = unique([yl(1) 0 yl(2)]);
      text(xlim*[1.02;-.02],ylim*[.98;.02], ...
           sprintf('%0.1f %s',h.YTick(1),y_unit), ...
           'VerticalAlignment','Bottom','HorizontalAlignment','right')
      text(xlim*[1.02;-.02],ylim*[.02;.98], ...
           sprintf('%0.1f %s',h.YTick(end),y_unit), ...
           'VerticalAlignment','Top','HorizontalAlignment','right')    
      h.YTickLabel = []; grid on
      set(gca,'UserData',kk)
    end

    %% Axis 2 - Coefficients
    axes('Position',p2+[0 (nK - kk)*p2(4) 0 0])
    imagesc(dat.x(1:nX), dat.ori(1:nX:end),reshape(dat.y_all(:,kk),nX,[])')
    colorbar
    
    if kk == nK, set(gca,'XTickLabel', strcat(get(gca,'XTickLabel'),' µm'))
    else         set(gca,'XTickLabel',{})
    end
    
    set(gca,'YTick',dat.ori(1:nX:end),'TickLength',[0 0])
    set(gca,'YTickLabel',strcat(get(gca,'YTickLabel'),'°'))
    set(gca,'UserData',kk)
    
    %% Axis 3 - Radon Transform
    axes('Position',p3+[0 (nK - kk)*p2(4) 0 0])
    
    if any(named('-al')), dat.algorithm = get_('-al');
    elseif any(named('alg')), dat.algorithm = get_('alg');
    end
    
    dat.y = dat.y_all(:,kk) - dat.y_base(kk);

    try dat = analysis.inverseRadon(dat);
    catch err, warning(err.getReport), continue
    end

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


if any(named('-latent'))
    %% Latent axes - explained variance 
    
    axes('Position',p3+[0 (nK - kk-1)*p2(4)+0.02 0 -0.035])
    if isfield(dat.latent), latent = dat.latent;
    else latent = get_('-latent');
    end
    explained = cumsum(latent)/sum(latent)*100;

    plot(1:numel(latent), explained, 'LineWidth',1.3)
    hold on, plot([nK nK nan xlim],[ylim nan explained(nK)*[1 1]],'LineWidth',1.2)
    xlim([1 numel(latent)]), title('\rmVariance Explained')
    set(gca,'YTickMode','manual','YTickLabel',strcat(get(gca,'YTickLabel'),'%'),'FontSize',7)
end

if nargout == 0, clear, end

end
