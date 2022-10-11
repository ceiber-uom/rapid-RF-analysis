
function totalSensitivity(dat, varargin )
% plots.totalSensitivity( data, ... )
% plots.totalSensitivity( data, time-points, ... )
% 
% Options: 
%  -image [radon data] : supply also the output of plots.plot_Radon_IMG
%                        (otherwise is recomputed)
%  -t [time-points] : alternate syntax
%  -row-size [4]    : number of subplots per row
%  -raw             : [todo] run radon reconstruction on mV value at
%                     timepoint (ignore PCA components) 
%  -interactive     : allow timepoint selection interactively by clicking 
%                     on waveform plot
% 
% v0.1 - 5 September 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

do_interactive = any(named('-in')); 

nK = size(dat.response_waves,2); 

if nargin > 1 && isnumeric(varargin{1}), timepoints = varargin{1};
    [~,timepoints] = arrayfun(@(t) min(abs(dat.time-t)), timepoints);        
elseif do_interactive,   [~,timepoints] = min(abs(dat.time)); 
elseif any(named('-t')), timepoints = get_('-t');
    [~,timepoints] = arrayfun(@(t) min(abs(dat.time-t)), timepoints);
else
    % selection based on maximum wave amplitude (in complex sense) of each
    % PCA (or NNMF, if so inclined) component

    % only consider maximal points in this window (could have edge effects)
    window = (dat.time > -0.1 & dat.time < 0.9*max(dat.time))'; 
    [~,timepoints] = arrayfun(@(k) max(abs(hilbert(dat.response_waves(:,k))) ... 
                                          .* window), 1:nK ); 
    [~,tzi] = min(abs(dat.time)); 
    timepoints = sort([tzi timepoints]); 
end

plots.standardFigure('Name','Total RF at timepoint'), clf

dat.response_baseline = mean(dat.response_waves(dat.time<=0, :));
dat.response_waves = dat.response_waves - dat.response_baseline;

% is_imp_s = (dat.time(1) == dat.time(2));
is_imp_s = isfield(dat,'psth');
if is_imp_s, y_unit = 'imp'; 
elseif any(named('-units-V')), y_unit = 'V';
else y_unit = 'mV';
  if max(abs(dat.response_waves(:))) < 0.5
    dat.response_waves = 1e3*dat.response_waves; 
    disp('Converting wave units to mV (assumed: V)')
  end
end

% Get radon transform data (or try to at any rate) 
if any(named('-im')), rdat = get_('-im');
elseif any(named('-raw')), rdat = []; 
else f = figure; rdat = plots.plot_radon_IMG(dat); delete(f);
end

%% Plot the component waves from which the timepoints were derived

npx = 4; 
if any(named('-row')), npx = get_('-row'); end
npy = ceil(numel(timepoints)/npx); 
sp_offset = 0; 

if isfield(dat,'response_waves')
    
    npy = npy + 1; 
    subplot(npy,1,1)
    sp_offset = npx; 

    for kk = 1:nK
        plot(dat.time, dat.response_waves(:,kk),'UserData',kk, ...
                                                'Hittest','off')
        hold on   
    end
    for ss = 1:dat.nStimuli % add stim bars to the waves plot
      rectangle('Position',dat.stim_bar(ss,0.1), ... 
                'FaceColor',[0 0 0 0.3], 'EdgeColor','none')
    end
    
    axis tight
    try tidyPlotForIllustrator, end %#ok<TRYNC>
    if do_interactive, hobj = gca; 
    else
        set(gca,'XTick',unique(round(dat.time(timepoints),2)), ... 
                'XTickLabelRotation',-90)
    end
    set(gca,'Position',get(gca,'Position') + [0 2 0 -1]/50)
    xlabel('time, s')
    if isfield(dat,'psth')
        ylabel('imp/s')
    else
        ylabel('mV')
    end
end

%% 

view_sinogram = ~isfield(rdat,'image') || any(named('-sino')); 

if do_interactive
    npx = 1; sp_offset = 1;
    plot(rdat.time(timepoints([1 1])), ylim,'Color',[0 0 0 0.3], ...
                                            'UserData','y-cursor');
    set(hobj,'ButtonDownFcn',@(a,b) on_timebase_click(a,b,dat,rdat));

    p = get(gca,'Position') ./ [1 1 1 nK]; 
    p(1) = p(1)+1.01*p(3);

    C = lines(max(7,nK));
    
    for kk = 1:nK
        uicontrol('style','checkbox','string',sprintf('PCA %d',kk), ...
                  'ForegroundColor',C(kk,:), 'Value',1, ...
                  'BackgroundColor','w', 'FontSize',11, ...
                  'units','normalized', 'UserData', kk, ...
                  'Position',p+[0 p(4)*(nK-kk) 0 0], ...
                  'Callback',@(a,~) on_timebase_click(a,[],dat,rdat))
    end
    
    % Scroll bar
    pos=get(hobj,'position');
    Newpos=[pos(1) pos(2)-0.12 pos(3) 0.03];
    h = uicontrol('style','slider','units','normalized','position',Newpos,...
    'Min',dat.time(1),'Max',dat.time(end),'SliderStep',[0.01,0.02],'callback',...
        @(a,~) S(a,[],dat,rdat));

end

for tt = 1:numel(timepoints) % show total RF at each timepoint

    total_rf_at_tt = determine_total_RF(dat, rdat, timepoints(tt), ...
                                                   view_sinogram);
    subplot(npy,npx,tt + sp_offset)
    if view_sinogram
      sino_axes = {unique(rdat.x), unique(rdat.ori)};
      imagesc(sino_axes{:}, total_rf_at_tt','UserData','sinogram')
      tidyPlotForIllustrator, axis tight xy
      set(gca,'YTick',sino_axes{2})
      set(gca,'YTickLabel',strcat(get(gca,'YTickLabel'),'°'))
    else 
      imagesc(rdat.range,rdat.range, total_rf_at_tt,'UserData','rf-map')
      axis image xy off
    end
    set(gca,'UserData',timepoints(tt))
    title(sprintf('t = %0.2f', dat.time(timepoints(tt))))
    set(gca,'Position',get(gca,'Position') + [-1 -1 2 2]/50)
end

h = flipud(findobj(gcf,'type','axes'));
h(cellfun(@isempty,{h.UserData})) = []; 
set(h,'CLim',[-1 1] * max(abs([h.CLim])))

c = interp1((-5:5)', redbluecmap, linspace(-5,5,101)); 
colormap(c), 

if npx > 1
    ch = colorbar('southoutside'); 
    ch.Position = [0.11 0.05 0.8 0.015];
    xlabel(ch,y_unit)
else
    ch = colorbar; ylabel(ch,y_unit)
end


% In principle you could compute the radon transform of the measured data 
% at these time-points in particular to compare to this visualisation,
% which is based repeated rounds of baseline subtraction. I'm leaving this
% as an exercise to the next person to tackle this (which might very well
% be me, some time in the future)

function S( src, ev, dat, r )

[~,timepoint] = min(abs(dat.time - src.Value));

h = findobj(gcf,'userdata','rf-map');
do_sinogram = isempty(h);

total_rf_at_tt = determine_total_RF(dat, r, timepoint, do_sinogram);

if do_sinogram
     h = findobj(gcf,'userdata','sinogram'); 
     set(h,'CData',total_rf_at_tt')
else set(h,'CData',total_rf_at_tt)
end

ok = isfinite(total_rf_at_tt); 
of_ = @(x) abs(reshape(x,[],1));

h.Parent.CLim = [-1 1] * max(of_(total_rf_at_tt(ok)), [],'omitnan');
title(h.Parent,sprintf('t = %0.2f', dat.time(timepoint)))
h.Parent.UserData = timepoint;

h = findobj(gcf,'userdata','y-cursor');
h.XData(:) = dat.time(timepoint);

return


%% UI interactivity functions  

function on_timebase_click( src, ev, dat, r )

if isempty(ev)

    ax = findobj(gcf,'type','axes');
    timepoint = [ax.UserData]; 
    timepoint = timepoint(1);

    ax = ax(cellfun(@isempty,{ax.UserData}));
    h = findobj(ax,'UserData',src.UserData);

    if src.Value, h.LineStyle = '-';
    else h.LineStyle = '--';
    end

else
    [~,timepoint] = min(abs(dat.time - ev.IntersectionPoint(1))); 
end


h = findobj(gcf,'userdata','rf-map');
do_sinogram = isempty(h);

total_rf_at_tt = determine_total_RF(dat, r, timepoint, do_sinogram);

if do_sinogram
     h = findobj(gcf,'userdata','sinogram'); 
     set(h,'CData',total_rf_at_tt')
else set(h,'CData',total_rf_at_tt)
end

ok = isfinite(total_rf_at_tt); 
of_ = @(x) abs(reshape(x,[],1));

h.Parent.CLim = [-1 1] * max(of_(total_rf_at_tt(ok)), [],'omitnan');
title(h.Parent,sprintf('t = %0.2f', dat.time(timepoint)))
h.Parent.UserData = timepoint;

h = findobj(gcf,'userdata','y-cursor');
h.XData(:) = dat.time(timepoint);

return


function on_image_click( hobj, ev, dat, r ) %#ok<DEFNU> 

error implement_this_feature

% ev.IntersectionPoint(1:2); 

% for each ori select bar closest to click

% validate that we're looking at the correct spot in terms of our
% selections

bar_selection = []; % passes 

u_ori = unique(r.ori(:)); 

color = hsv(numel(u_ori)); 
color = color ./ sqrt(sum(color,2)) .* [.9 .85 1];
color(:,4) = 0.5; % alpha 

plots.standardFigure('Name','Responses to selected stimuli'), clf

for pp = bar_selection

    ori_id = (u_ori == this_ori); 
    plot(dat.time, dat.raw_waveform(pp,:),'Color',color(ori_id,:))

end

tidyPlotForIllustrator, xlim(dat.time([1 end]))


function img = determine_total_RF(dat, rdat, idx, sinogram)

h = findobj(gcf,'Style','checkbox');
nK = size(dat.response_waves,2); 
if isempty(h), ok = true(1,nK);
else 
    [~,seq] = sort([h.UserData]); 
    ok = [h(seq).Value] ~= 0;
end

if sinogram

  nO = numel(unique(rdat.ori)); 
  sinogram_ = @(y) reshape(y,[],nO);
  img = sinogram_(rdat.y_all(:,ok) * dat.response_waves(idx,ok)');

elseif isempty(rdat)
  error TODO_implement_RAW_mode    
else
  img = zeros(size(rdat.images{1})); 
  for k = find(ok)
    img = img + rdat.images{k} * dat.response_waves(idx, k); 
  end
end

