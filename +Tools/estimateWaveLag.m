function timing = estimateWaveLag(wave, varargin)
% timing = estimateWaveLag(wave, [time, expoData], [...] )
%   estimate the response latency of the input wave to a step stimulus
%   using a lagged mutual-information-type estimator. 
% 
% Optional inputs:
%   -plot, -tf ### (Defaut: from ExpoData), -roi ### (Default: 0.5 s)

named = @(n) strncmpi(varargin,n,length(n));
do_plot = any(named('-plot')); 

%%

if nargin == 0 || nargin > 1e3, do_plot = true; 
    % nargin > 1e3 == script-mode
     
     disp('Running demonstration (PCA 1):')
    
     Tools.loadPhysiology('8.6.1 #13','-plot','-pca','-nK',6);
     wave = score(:,1); %#ok<NODEF>     
     
     % Remove variables that non-test-case use doesn't have access to
     clear coeff coe_scale filename hekaData nPasses nStimuli 
     clear passes resting score stim_bar
     
elseif nargin < 3
     expoData = evalin('caller','expoData');
     time = evalin('caller','time');
else
     expoData = varargin{2};
     time = varargin{1};
end

if exist('named','var') && any(named('-tf')) % Usage: -tf 2 (Defaut: from ExpoData)
     tf = varargin{find(named('-tf'))+1}; 
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

if exist('named','var') && any(named('-roi')) % Usage: -roi 0.4 (default)
     lags = varargin{find(named('-roi'))+1}; 
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

if do_plot
    plot(lags, 10 * MI,'LineWidth',1.2,'Color',C(2,:))
    plot(time(t_MI),wave(t_MI),'*','Color',C(1,:),'LineWidth',1.3)
    plot(time([1 t_MI]), [1 1]*y_ref, '--','Color',C(3,:),'LineWidth',1.1)    
end

t_os = find(wave(1:t_MI)*sign(wave(t_MI)) < y_ref*sign(wave(t_MI)),1,'last');
t_os = round(max(t0, [t_os t_MI]*[1.5;-0.5])); 

if do_plot
    plot([1 1]*time(t_os),ylim,'-','Color',C(3,:),'LineWidth',1.1)
end

%% Export to "timing" structure for output

timing.zero = t0; 
timing.index = t_os; 
timing.latency = time(t_MI); 
timing.duration = sum(input == 2); % number of time samples for which the stimulus is on (1/2 cycle)

timing.nStimuli = round(max(time) * tf);

timing.MI_lag = lags;
timing.MI_val = MI;


%% Direct plot of membrane potential    



% clear t0 t tf dt in_lagged mi MI input lags nT ok p_i p_ij ss t_MI t_os tt y_bins y_ref
