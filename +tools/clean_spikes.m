
function clean_spikes(data, varargin)
% tools.clean_spikes(data, [pass_ids], ... )
% 
% The spike detection can sometimes muck up for spike trains where Vm
% doesn't return down to the resting membrane potential (or, more to the
% point, below the spike detection threshold) between spikes. 
% 
% Credit to EB for spotting this bug. 
% 
% This fixes this by throwing out any spikes whose waveforms aren't
% approximately (90%) max at t=0 (about sample 16)
% 
% Options: 
% -q [0.95] : the sample at t=0 must be greater than this fraction of samples
% -z [16]   : which sample is t=0?
% -p [ids]  : only process which pass_ids? 
% -no-p     : do not plot
% -no-w     : do not recompute PSTH
% 
% v0.1 - Calvin Eiber - 8 July 2023

%% Get inputs

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

do_plot = ~any(named('-no-p'));

pass_ids = 1:data.nPasses; 
if any(named('-p')), pass_ids = get_('-p');
elseif nargin > 1 && isnumeric(varargin{1}), pass_ids = varargin{1};
end

zero_sample = 16;
if any(named('-z')), zero_sample = get_('-z'); end

min_quantile = 0.95; 
if any(named('-q')), min_quantile = get_('-q'); end


%%

spk = data.hekaData.Spikes;

time = data.wave_time;
passes = ismember(spk.passIdx,pass_ids);

val = mean(spk.spikeWaveforms <= spk.spikeWaveforms(zero_sample,:),1);
cull = (val < min_quantile)  & passes;

if do_plot
    %%
    
    clf
    subplot(2,2,1)
    plot( spk.spikeWaveforms(:, passes & ~cull ), '-', 'color', [0 0 0 0.1])
    try tidyPlotForIllustrator, end %#ok<TRYNC> 
    title('Good spikes')
    
    subplot(2,2,2)
    plot( spk.spikeWaveforms(:, passes & cull ), '-', 'color', [0.5 0 0 0.1])
    try tidyPlotForIllustrator, end %#ok<TRYNC> 
    title('Bad spikes')
    
    subplot(2,2,[3 4])
    plot( time(spk.timeIdx(passes & ~cull)), ...
               spk.passIdx(passes & ~cull),'k.'), hold on
    plot( time(spk.timeIdx(passes & cull)), ...
               spk.passIdx(passes & cull),'rx')
    
    axis tight, xlim(data.wave_time([1 end]))
    try tidyPlotForIllustrator, end %#ok<TRYNC> 
end

%% Remove bad spikes from data

data.hekaData.Spikes.timeIdx(cull) = []; 
data.hekaData.Spikes.passIdx(cull) = []; 
data.hekaData.Spikes.timeSeconds(cull) = []; 
data.hekaData.Spikes.spikeWaveforms(:, cull) = []; 

%% Recompute PSTH
if any(named('-no-w')), return, end

bins = arrayfun(@(t) find(time == t), data.time); 
wave = data.psth.wave * median(diff(data.psth.time)); 

for ii = reshape(pass_ids,1,[]) % make PSTH 
   pass = data.hekaData.Spikes.passIdx == ii;
    if ~any(pass), continue, end
    wave(:,ii) = hist(data.hekaData.Spikes.timeIdx(pass),bins); %#ok<HIST> 
end

data.psth.wave = wave / median(diff(data.psth.time));

%%
return


