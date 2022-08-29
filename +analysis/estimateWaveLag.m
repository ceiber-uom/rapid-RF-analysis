function timing = estimateWaveLag(wave, varargin)
% timing = estimateWaveLag(wave, [time, expoData], [...] )
%   estimate the response latency of the input wave to a step stimulus
%   using a lagged mutual-information-type estimator followed by an
%   examination of departure from the 95% PI of the first-guess baseline
%   region
% 
% Optional inputs:
%   -plot, -tf ### (Defaut: from ExpoData), -roi ### (Default: 0.5 s)
% 
% V0.2 - 29 August 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

do_plot = any(named('-plot')); 

%%

if nargin == 0 || nargin > 1e3, do_plot = true; 
    % nargin > 1e3 == script-mode
     
     disp('Running demonstration (PCA 1):')
    
     utils.load('8.6.1 #13','-plot','-pca','-nK',6);
     wave = response_waves(:,1) * 1e3; %#ok<NODEF>     
     
     % Remove variables that non-test-case use doesn't have access to
     clear activations response_waves response_scaleFactor filename
     clear hekaData nPasses nStimuli passes stim_bar resting_potential
elseif isstruct(wave) 
    timing = default_dataset_script(wave, varargin{:});
    return

elseif nargin < 3


     expoData = evalin('caller','expoData');
     time = evalin('caller','time');
else
     expoData = varargin{2};
     time = varargin{1};
end

if any(named('-tf')) % Usage: -tf 2 (Defaut: from ExpoData)
     tf = get_('-tf'); 
elseif isstruct(expoData), 
     tf = expoData.passes.events{2}.Data{1}(3); 
else tf = expoData; clear expoData
end

% Construct input waveform
input = 0*time; 

% n_ii = nStimuli; % Compute total latency for entire trace (unstable!)
% n_ii = 1;        % Compute latency just for first rising and falling edge

t = [0 0.5]/tf; 
input(time > t*[1;0]) = 2;
input(time > t*[1;1]) = 1;
input(time > t*[1;2]) = -1;

nT = sum(input >= 0); 
t0 = find(time > 0, 1); 

if do_plot
    C = lines(7); G = @(v) [v v v]/10;    
    clf, plot(time(input>=0),input(input>=0),'-','LineWidth',1.7,'Color',G(3))
    hold on, plot(time,wave,'Color',C(1,:)), pause(0.01)
    tidyPlotForIllustrator,  xlim(time([1 end]))
end

if any(named('-roi')) 
     lags = get_('-roi'); 
else lags = 0.4; % maximum latency
end

y_bins = conv(linspace(min(wave),max(wave),65),[1 1]/2,'valid'); 
lags = 0:round(lags / median(diff(time))); 
MI = zeros(size(lags)); 

%% Calcuate lagged mutual information
for tt = 1:length(lags) % compute mutual information between input, lagged response 
    
    dt = lags(tt);     
    in_lagged = [zeros(1,dt) input(1:end-dt)];
        
    if in_lagged(end) < 0       
        ok = (1:nT)+dt;
        p_i = hist(wave(ok),y_bins) / nT;

        for ss = 0:max(input)
            % snip = ; 
            p_ij = hist(wave(ok(input(1:nT) == ss)),y_bins) / nT;
            mi = nansum(p_ij .* log(p_ij./p_i./sum(input == ss)*nT));
            MI(tt) = MI(tt) + mi;
        end
    else
        if ~isempty(ok), 
            ok = (1:nT)+dt;
            p_i = hist(wave(ok),y_bins) / nT;
            ok = [];
        end
        in_lagged(1:(end-nT)) = -2; % Trim off leading edge
        
        for ss = 0:max(input)
            % snip = score(in_lagged == ss,kk);
            p_ij = hist(wave(in_lagged == ss),y_bins) / nT;
            mi = nansum(p_ij .* log(p_ij./p_i./sum(input == ss)*nT));        
            MI(tt) = MI(tt) + mi;
        end

%         MI(tt:end) = []; lags(tt:end) = []; 
%         break
    end
end


%% Get latency (and beginning of trace) from MI maximum
lags = lags * mean(diff(time)); 

[~,t_MI] = max(MI);
[~,t_MI] = min(abs(time - lags(t_MI)));
y_ref = median(wave(t0:t_MI));

pred_interval = []; 


if any(named('-do-pi')) % no predint correction
    %% Looks like the latency is somewhat overestimated by MI

    predint_roi = 0.75;     
    predint_alpha = 0.05; 

    if any(named('-pi-r')), predint_roi   = get_('-pi-r'); end
    if any(named('-pi-a')), predint_alpha = get_('-pi-a'); end

    predint_q = [predint_alpha/2 1-predint_alpha/2]; % quantiles 

    predint_roi = ceil(predint_roi*t_MI); 

    pred_interval = quantile(wave(1:predint_roi), predint_q);    
    in_ci = (wave >= pred_interval(1) & wave <= pred_interval(2)); 

    if do_plot

        fill(time([1 t_MI t_MI 1 1]), pred_interval([1 1 2 2 1]), ...
             C(6,:),'FaceAlpha',0.2,'EdgeColor','none')
        plot(time(t_MI),wave(t_MI),'o','Color',C(1,:), ...
                                   'LineWidth',1,'MarkerFaceColor','w')
    end

    cix = find( ~in_ci(1:t_MI) & [true; in_ci(1:t_MI-1)]); 
    if any(cix > predint_roi), t_MI = cix(end); end

end

if do_plot
    plot(lags, 10 * MI,'LineWidth',1.2,'Color',C(3,:))
    plot([1 1]*time(t_MI),ylim,'-','Color',C(3,:) ) 
%     plot(time(t_MI),wave(t_MI),'v','Color',C(1,:),'LineWidth',1.3,...
%                                             'MarkerFaceColor',C(1,:))
    plot(time([1 t_MI]), [1 1]*y_ref, '-','Color',C(3,:),'LineWidth',1.1)    
end

t_os = find(wave(1:t_MI)*sign(wave(t_MI)) < y_ref*sign(wave(t_MI)),1,'last');
% t_os = round(max(t0, [t_os t_MI]*[1.5;-0.5])); 

if do_plot
    plot(time(t_os),wave(t_os),'o','Color',C(1,:),'LineWidth',1.1,'MarkerFaceColor',C(1,:))
end

%%

%% Export to "timing" structure for output


timing.index = t_os; 
timing.latency = time(t_os); 
timing.duration = sum(input == 2); % number of time samples for which the stimulus is on (1/2 cycle)

timing.zero_index = t0; 
timing.MI_index = t_MI;
timing.MI_latency = time(t_MI);

timing.nStimuli = round(max(time) * tf);

timing.MI_lag = lags;
timing.MI_val = MI;

if ~isempty(pred_interval)
    timing.PI_val = pred_interval;
    timing.PI_roi = [1 predint_roi];
    timing.PI_alpha = predint_alpha;
end

return


function out = default_dataset_script(dat, varargin)

named = @(n) strncmpi(varargin,n,length(n));


out = []; 
for yy = 1:size(dat.response_waves,2)

    if any(named('-plot')), figure(yy), clf, end

    t = analysis.estimateWaveLag(dat.response_waves(:,yy) * 1e3, ...
                                 dat.time, ...
                                 dat.expoData, ...
                                 varargin{:});

    if isempty(out), out = t;
    else out(end+1) = t;
    end
end