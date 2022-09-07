
function totalSensitivity(dat, varargin )
% plots.totalSensitivity( data, ... )
% plots.totalSensitivity( data, time-points, ... )
% 
% Options: 
%  -image [radon data] : supply also the output of plots.plot_Radon_IMG
%                        (otherwise is recomputed)
%  -t [time-points] : alternate syntax
%  -row-size [4] : number of subplots per row
%  -raw 
% 
% v0.1 - 5 September 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

nK = size(dat.response_waves,2); 

if nargin > 1 && isnumeric(varargin{1}), timepoints = varargin{1};
    [~,timepoints] = arrayfun(@(t) min(abs(dat.time-t)), timepoints);        
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

%% Plot the component waves from which the timepoints were derived

npx = 4; 
if any(named('-row')), npx = get_('-row'); end
npy = ceil(numel(timepoints)/npx); 
sp_offset = 0; 

if isfield(dat,'response_waves')
    
    npy = npy + 1; 
    subplot(npy,1,1)
    sp_offset = npx; 

    for k = 1:nK
        plot(dat.time, 1e3*dat.response_waves(:,k)), hold on   
    end
    for ss = 1:dat.nStimuli % add stim bars to the waves plot
      rectangle('Position',dat.stim_bar(ss,0.1), ... 
                'FaceColor',[0 0 0 0.3], 'EdgeColor','none')
    end
    
    axis tight
    try tidyPlotForIllustrator, end %#ok<TRYNC> 
    set(gca,'XTick',unique(round(dat.time(timepoints),2)),'XTickLabelRotation',-90)
    set(gca,'Position',get(gca,'Position') + [0 2 0 -1]/50)
end

%% 

if any(named('-im')), rdat = get_('-im');
elseif any(named('-raw')), rdat = []; 
else 
    f = figure; rdat = plots.plot_radon_IMG(dat); delete(f);
end


for tt = 1:numel(timepoints) % show total RF at each timepoint

    total_rf_at_tt = 0 * rdat.image; 

    for k = 1:nK 
        total_rf_at_tt = total_rf_at_tt + rdat.images{k} * ...
                        dat.response_waves(timepoints(tt), k); 
    end

    subplot(npy,npx,tt + sp_offset)
    imagesc(rdat.range,rdat.range, 1e3*total_rf_at_tt)
    axis image xy off
    title(sprintf('t = %0.2f', dat.time(timepoints(tt))))
    set(gca,'Position',get(gca,'Position') + [-1 -1 2 2]/50)
end

h = flipud(get(gcf,'Children')); h(1) = []; 
set(h,'CLim',[-1 1] * max(abs([h.CLim])))

c = interp1((-5:5)', redbluecmap, linspace(-5,5,101)); 
colormap(c), 

ch = colorbar('southoutside'); 
ch.Position = [0.11 0.05 0.8 0.015];
xlabel(ch,'mV')

% In principle you could compute the radon transform of the measured data 
% at these time-points in particular to compare to this visualisation,
% which is based repeated rounds of baseline subtraction. I'm leaving this
% as an exercise to the next person to tackle this (which might very well
% be me, some time in the future)

