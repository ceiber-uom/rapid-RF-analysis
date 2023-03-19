
function data = load(varargin)
% tools.load( filename, ... )
% 
% Select a .mat file from ./MAT/ containing hekaData and expoData
%   (as generated by convert_HEKA_to_MAT), load the selected file, and
%   perform some standard analyses. 
% 
% filename may be a full path to a file, a filename in ./MAT, or an
%   abbreviated file specification e.g. "8.6.1 #9" or "20180806_Cell_1 #9"
% 
% The default data directory is ./MAT; this can be changed by manually
%  selecting a different directory or by setting -dir '/path/to/new/dir'.
% Note that the default data directory persists between calls to this
%  function, so for repeatable analysis scripts the following example
%  usage is recommended: 
% 
% ` list = dir('path/to/my/*.mat) % list of files to analyse
% ` p_ = @(x) [x(1).folder filesep x(1).name]; % path expander
% ` 
% ` for ii = 1:numel(list) 
% `     data = tools.load(p_(list(ii)), options )
% `     do_some_analysis( data, ... )
% ` end
% 
% Load & pre-process options: 
% -dir './MAT' : set data directory
% -spike : remove action potential waveforms from the membrane potential
%          using analysis.removeActionPotentials (enabled by default)
% -psth  : compute PSTH of recorded spikes. Related options: 
%          -bin [100] : Set bin size (in samples)
% -f01   : compute zeroth- and first-harmonic response to periodic stimulus.
%           Not implemented in conjunction with -psth, use expo tools. 
% -pca   : compute PCA of waves using analysis.linear, see below.
% -nnmf  : compute NNMF of spike rate responses using analysis.linear
% -dump  : instead of returning a data structure, dump the data into the
%          calling workspace (default if nargout = 0)
% 
% Response factorisation options: 'pca', 'nnmf', 'ica'. These correspond to
%   the working modes of analysis.linear, see doc analysis.linear. 
% 
% If 'X' is a is a (time-by-stimuli) array of membrane potentials or 
%   spike-rates, a linear decomposition of 'DATA' is a pair of matrices 
%  'activations' (stimuli-by-nK) and 'response_waves' (time-by-nK) which
%   satisfies the equation: 
% 
%          X = (Y.response_waves * Y.activations') + Y.baseline
% 
% This analysis splits the input data into components which capture the
%   different features of the overall response. for PCA, these are the
%   directions of maximum variance. The following options are passed to
%   analysis.linear: 
%    -nK [6] : set default number of components to be returned. 
%    -rest   : disable baseline subtraction (not recommended) 
% 
% Basic visualisation: 
% -plot : show loaded membrane potential traces and 
% -plot-spectra : debug plot of spectra of responses (requires -f01)
% -plot-f1      : debug scatterplot of max(ifft) vs f1 (requires -f01)
% 
% Persistent analysis options can be saved in the .mat file as fields in
% the .options structure to customise the default analysis of each file. 
% 
% options.apply_offset [offset in ms] - if this exists, shift all response
%    waves by [options.apply_offset] ms (corrects for extremely long
%    latency responses bleeding into next trail) 
% options.apply_bin_size [bin size in µs] - modify the default value of
%    -bin for this file (default 100, can decrease if spikerates are high
%    and lots of spikes per bin or, more likely, increase if few/no spikes
%    observed). 
% 
% Output: 
% 
% .filename - input filename with underscores escaped ('_' => '\_')
% .hekaData - HEKA data for input file
% .expoData - EXPO data for input file
% .time     - peri-stimulus time in seconds
% .passes   - vector of block IDs for each stimulus. 
% .nPasses  - total number of passes (nP) 
% .nStimuli - how many stimulus presentations were there per pass?
%             (related to temporal frequency of stimulus)
% .stim_bar - utility function for showing stimulus in context with trace.
% 
% [if -psth requested]: 
% .psth.time - time vector for PSTH waveform (length = nB)
% .psth.wave - [nB x nP] matrix of per-trial binned spike-rates 
% 
% [if -f01 requested]: 
% .f0      - [real] time-averaged post-stimulus membrane potential
% .f1      - [complex] f1 Vm response to stimulus at stim frequency
% .f1_wave - sine wave showing f1 response (in V)
% .f1_roi  - identified post-stimulus region of time
% 
% [if -pca or -nnmf requested]: 
% .activations: [nP x nK] value of each component for each stimulus. 
% .response_waves: [nT × nK] response waveform for each component 
% .resting_potential: pre-stimulus potential 
% .response_scaleFactor: scale factor to achieve unit activation scores. 
%                        (optional, if -get-s set)
% 
% Example usage of data.stim_bar: 
% ` plot(data.time, data.hekaData.PassData, 'Color', [0 0 0 0.05])
% ` xlim([min(data.time) max(data.time)])
% ` 
% ` for ss = 1:data.nStimuli,
% `   rectangle('Position',data.stim_bar(ss,0.1), ... 
% `             'FaceColor',[0 0 0 0.5], 'EdgeColor','none')
% ` end
% 
% Version 2 - 28 August 2022 - Calvin Eiber <ceiber@ieee.org>
%                              Refactored from old Tools.loadPhysiology

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

filename = parse_filename( varargin ); 

disp(['Loading ' filename]), 
load(filename,'hekaData','expoData','options');

%% Get basic information about stim and recorded data 

fs = 1./hekaData.SweepHeader.SampleInterval;
passes = double(expoData.passes.BlockIDs(2:2:end));
time = (0:size(hekaData.PassData,1)-1) / fs ;

delay = diff(expoData.passes.StartTimes);
delay = mean(delay(1:2:end)) / 1e4; 
time = time - delay; 

duration = double(expoData.passes.EndTimes - ...
                  expoData.passes.StartTimes) / 1e4;
duration = mean(duration(2:2:end));

%% Remove Action Potentials from membrane current
if ~any(named('-spike'))
    hekaData = analysis.removeActionPotentials(hekaData); 
end

%% If needed, apply a time offset (if response has crept into next pass)
if isfield(options,'apply_offset') && ~isempty(options.apply_offset)    
    [hekaData,expoData,time] = do_apply_offset(hekaData,...
                                               expoData,time,options);
end

%% Direct plot of membrane potential    

stim_rate = expoData.passes.events{2}.Data{1}(3);
stim_bar = @(i,h) [(i-1)/stim_rate ylim*[h;1-h] 1/stim_rate/2 ylim*[-1;1]*abs(h/3)];
nStimuli = round(duration * stim_rate);  

if any(named('-plot'))
    
    figure(1), clf
    set(gcf,'Color','w','Position',[50 560 960 330],'Name','Membrane Potential')

    plot(time, 1000*hekaData.rawPassData, 'Color', [0 0 0 0.05])
    xlim([min(time) max(time)])

    for ss = 1:nStimuli
        rectangle('Position',stim_bar(ss,0.1),'FaceColor',[0 0 0 0.5], 'EdgeColor','none')
        % rectangle('Position',stim_bar(ss+0.5),'FaceColor',[0  0  0  .4], 'EdgeColor','none')
    end

    set(gca,'XTickMode','manual','XTickLabel',strcat(get(gca,'XTickLabel'),' s'))
    set(gca,'YTickMode','manual','YTickLabel',strcat(get(gca,'YTickLabel'),' mV'))
    set(gca,'TickDir','out','box','off','Color','none')

    if ~isempty(which('tidyPlotForIllustrator')), tidyPlotForIllustrator, end
end

%% Begin constructing output data structure

[~,fn,~] = fileparts(filename); 

data = struct; 
data.filename = strrep(fn,'_','\_');
data.hekaData = hekaData;
data.expoData = expoData;

data.time = time; 
data.passes = passes;
data.nPasses = numel(passes);
data.nStimuli = nStimuli;
data.stim_bar = stim_bar; % utility finction 

%% Compute PSTH (if requested)
if any(named('-psth'))

    wave_time = time; % Possible "jumps" - ignore for analysis 
    
    if any(named('-bin')), time = get_('-bin'); 
    elseif isfield(options,'apply_bin_size')
         time = options.apply_bin_size; 
         fprintf('Using non-standard bin size %d (default: 100) for this file\n', time)        
    else time = 100;
    end
    
    time = 0:time:size(hekaData.PassData,1);
    time = time(1:end-1)/2 + time(2:end)/2;

    wave = zeros(numel(time), size(hekaData.PassData,2));
    
     % 24.02.22 Increase spike detection threshold for cells that have a 
    % resting Vm or depolarisation exceeding -40 mV
    file = extractBefore(expoData.FileName,'.');
    threshold_15 = {'20210708_Cell_05#9[Radon_Flicker_ACH]';'20210708_Cell_05#14[Radon_Flicker_ACH]'};
    threshold_20 = {'20190904_Cell_02#15[Radon_Flicker_ACH]';'20190904_Cell_02#16[Radon_Flicker_ACH]';...
        '20190904_Cell_02#17[Radon_Flicker_ACH]';'20190904_Cell_02#18[Radon_Flicker_ACH]';...
        '20190904_Cell_02#2[off_size_ms]';'20210629_Cell_02#8[off_size_ms]';...
        '20210629_Cell_02#10[Radon_Flicker_ACH]';'20210629_Cell_02#13[Radon_Flicker_ACH]';...
        '20220228_Cell_02#9[Radon_Flicker_ACH]';'20220401_Cell_02#8[Radon_Flicker_ACH]';...
        '20220811_Cell_02#6[Radon_Flicker_ACH]';'20220929_Cell_01#6[Radon_Flicker_ACH_long]';...
        '20220929_Cell_01#5[on_size_ms]'}; 
    threshold_25 = {'20191023_Cell_05#7[Radon_Flicker_ACH]';'20191023_Cell_05#8[Radon_Flicker_ACH]';...
        '20210708_Cell_02#10[Radon_Flicker_ACH]';'20220923_Cell_01#7[Radon_Flicker_ACH_long]';...;...
        '20220923_Cell_01#6[on_size_ms]'}; 
    threshold_30 = {'20190904_Cell_02#12[sf_ms]';'20210708_Cell_02#12[Radon_Flicker_ACH]';...
        '20210708_Cell_02#11[on_size_ms]';'20220127_Cell_03#12[Radon_Flicker_ACH]';...
        '20220127_Cell_03#13[Radon_Flicker_ACH]';'20220228_Cell_02#7[on_size_ms]';'20220228_Cell_02#10[sf_ms]';...
        '20220805_Cell_02#6[Radon_Flicker_ACH]';'20220811_Cell_02#9[Radon_Flicker_ACH]';...
        '20220811_Cell_03#6[Radon_Flicker_ACH]';'20220811_Cell_03#8[Radon_Flicker_ACH]';...
        '20220812_Cell_01#6[Radon_Flicker_ACH_long]';'20220929_Cell_03#6[Radon_Flicker_ACH_long]'};
    threshold_35 = {'20220218_Cell_02#8[Radon_Flicker_ACH]'};
    threshold_45 = {'20200116_Cell_02#15[Radon_Flicker_ACH]'};
    
    if any(contains(threshold_15,file))
        idx = max(hekaData.Spikes.spikeWaveforms >= -0.015);
    elseif any(contains(threshold_20,file))
        idx = max(hekaData.Spikes.spikeWaveforms >= -0.02);
    elseif any(contains(threshold_25,file))
        idx = max(hekaData.Spikes.spikeWaveforms >= -0.025);
    elseif any(contains(threshold_30,file))
        idx = max(hekaData.Spikes.spikeWaveforms >= -0.03);
    elseif any(contains(threshold_35,file))
        idx = max(hekaData.Spikes.spikeWaveforms >= -0.035);
    elseif any(contains(threshold_45,file))
        idx = max(hekaData.Spikes.spikeWaveforms >= -0.045);
    end
    if any(contains(threshold_15,file)) || any(contains(threshold_20,file)) ...
       || any(contains(threshold_30,file)) || any(contains(threshold_25,file)) || ...
        any(contains(threshold_35,file)) || any(contains(threshold_45,file)) 
        hekaData.Spikes.spikeWaveforms = hekaData.Spikes.spikeWaveforms(:,idx); 
        hekaData.Spikes.timeSeconds = hekaData.Spikes.timeSeconds(:,idx);
        hekaData.Spikes.passIdx = hekaData.Spikes.passIdx(:,idx); 
        hekaData.Spikes.timeIdx = hekaData.Spikes.timeIdx(:,idx); 
    end


    for ii = 1:size(hekaData.PassData,1) % make PSTH 
        pass = hekaData.Spikes.passIdx == ii;
        if ~any(pass), continue
        else wave(:,ii) = hist(hekaData.Spikes.timeIdx(pass),time); %#ok<HIST> 
        end
    end

    time = wave_time(time); % If there's a jump, need to know
    % time = (time * hekaData.SweepHeader.SampleInterval) - delay;    

    data.psth.time = time;
    data.psth.wave = wave / median(diff(time));

    % EB 3-point smooth: makes RF map less sharp
    % No smooth for: 20220401_Cell_02, 20220513_Cell_02 (nps=1)
    if any(named('-smooth'))
        nps = 3; 
        pad = (nps-1)/2;
        smooth = @(y) conv2(y(:,[ones(1,ceil(pad)) 1:end end*ones(1,floor(pad))]), ...
             ones(1,nps)/nps,'valid');
        data.psth.wave = smooth(data.psth.wave);
    end
      
end

%% "Traditional" F0/F1 analysis of membrane potentials

if any(named('-f01'))
    
    if any(named('-psth'))
        error TODO_spike_f01
    end
    
    roi = [max(time)-min(1,2/stim_rate) inf]; 
    roi = (time >= min(roi) & time <= max(roi)); 

    % stim_trace = time > 0 & mod(time*stim_rate,1) < 0.5; 
    stim_trace = (time > 0) .* sin(2*pi*time*stim_rate); 

    S = fft(stim_trace(roi));  
    F = fft(hekaData.PassData(roi,:)) / sum(roi);
    S(1) = 0; % S = S ./ sqrt(sum(abs(S).^2));  
    F(1,:) = 0; 

    ok = 1:length(S)/2;

    f1_wave = ifft(((diag(S / sum(abs(S)))*F)),'symmetric')*sum(roi)*2; 
    f1 = (4*S(ok)*F(ok,:)) / sum(abs(S)); %  * sum(roi) * 2);

    if any(named('-plot-s'))
        %% Debug: spectrum of the stimulus S
        clf
        hz = linspace(0,1,sum(roi))*fs; 
         
        loglog(hz(ok),abs(S(ok)))
        axis tight, hold on
        if ~exist('pp','var'), [~,pp] = max(abs(f1)); end
        loglog(hz(ok),abs(F(ok,pp)))
        loglog(hz(ok),abs(F(ok,pp) .* S(ok)'))
        plot(hz(ok),abs(f1(pp))*ones(size(ok)),'--')
    end
    if any(named('-plot-f'))
        %% Debug: plot max of ifft wave vs computed f1
        clf
        plot(max(f1_wave),abs(f1),'.')
        hold on, plot(xlim,xlim,'Color',[0 0 0 0.3])
    end

    di = mean(exp(1i*angle(f1))); 
    f1 = f1 / (di/abs(di)) * -sign(real(di)); 
    f0 = mean(hekaData.PassData(roi,:));

    data.f0 = f0;
    data.f1 = f1;
    data.f1_wave = f1_wave;
    data.f1_roi = roi; 
end

%% Final standard analysis: wave decomposition (if requested)

if any(cellfun(@(x) any(named(x)), ... % any of these:
              {'pca','-pca','nnmf','-nnmf','ica','-ica'}))
   
    R = analysis.linear(data, varargin{:}); 

    data.activations = R.activations;
    data.response_waves = R.response_waves;
    data.resting_potential = R.baseline;
    
    if any(named('-psth')) && ~any(named('-exact-time'))
        % make sure that the time vector matches data.response_waves
        data.wave_time = data.time;
        data.time = data.psth.time;
    end
    
    if any(named('-get-s')) % optional, quality-of-life upgrade to omit
        data.response_scaleFactor = R.response_scaleFactor;
    end
end


if nargout == 0 || any(named('-dump'))
  %% Dump the loaded data into the caller workspace
  for f = fieldnames(data)', assignin('caller',f{1}, data.(f{1})); end
  clear data  
end

return
    
    
    


function filename = parse_filename( arg_in )
% from input arguments, determine the filename

named = @(n) strncmpi(arg_in,n,length(n));

persistent fp % file path
if isempty(fp) || any(fp == 0), 
    % default filepath: ./MAT
    % (relative to the directory contining +tools.load)
    fp = regexprep(fileparts(mfilename('fullpath')),'([/\\])\+.*','$1MAT$1');
    if ~exist('fp','dir')
      fp = regexprep(fileparts(mfilename('fullpath')),'([/\\])\+.*','$1data$1');
    end
end

if any(named('-dir')), fp = arg_in{find(named('-dir'),1) + 1}; end
if fp(end) ~= filesep, fp = [fp filesep]; end
if ~exist('fp','dir')

end

fn = cellfun(@(v) ischar(v) && any(v=='#'), arg_in); 

if any(fn), fn = arg_in{find(fn,1)};
  if exist(fn,'file'), filename = fn; return 
  elseif exist([fp fn],'file'), filename = [fp fn]; return
  else
    % Parse code of the form "8.6.1 #9" 
    % also can accept "20180806_Cell_1 #9"
    dc_code = str2double(regexp(fn,'\d+','match'));
    if ~exist(fp,'dir'), error('folder %s not found.', fp), end

    if numel(dc_code) == 3 && dc_code(1) > 100        
        dc_code = [mod(dc_code(1),[1e4 1e2]) dc_code(2:end)];
        dc_code(1) = (dc_code(1) - dc_code(2))/100;
    elseif numel(dc_code) ~= 4, error('%s not a valid code string.', fn), 
    end

    fn = sprintf('%02d%02d_Cell_%02d',dc_code(1:3));
    list = dir([fp '*'  fn '*.mat']);

    if isempty(list) % try just %d for cell_id
        fn = sprintf('%02d%02d_Cell_%d',dc_code(1:3));
        list = dir([fp '*'  fn '*.mat']);
    end
    
    if isempty(list), error('%s not found in %s', fn, fp), end
    hash_no = cellfun(@(s) str2double(regexp(s,'(?<=#\w*)\d+','match')), {list.name}); 
    
    if ~any(hash_no == dc_code(end)), 
        error('%s #%d not found in [%s%s].', fn, dc_code(end), sprintf('#%d,',hash_no),8)
    end

    fn = list(hash_no == dc_code(end)).name;
    filename = [fp fn]; return
  end
else
    [fn,fp] = uigetfile('*.mat','Select data',fp);    
    filename = [fp fn]; return    
end



function [hekaData,expoData,time] = do_apply_offset(hekaData,...
                                                    expoData,time,options)

% 1e4 is the (fixed) expo sample rate in Hz 

sel = time < options.apply_offset;
spk_select = sel(hekaData.Spikes.timeIdx); 
pass_delay = diff(expoData.passes.StartTimes(1:2:end));
pass_delay = mean(pass_delay) / 1e4;
   
spk_passID = hekaData.Spikes.passIdx; 
spk_passID(spk_select) = spk_passID(spk_select) - 1;
time(sel) = time(sel) + pass_delay; 
[time, order] = sort(time);

hekaData.PassData = hekaData.PassData(order,:); 
if isfield(hekaData,'rawPassData')
  hekaData.rawPassData = hekaData.rawPassData(order,:); 
end

% Invert order mapping to re-assign spike-times
order(order) = 1:length(order);

hekaData.Spikes.timeIdx = order(hekaData.Spikes.timeIdx); 
hekaData.Spikes.passIdx = spk_passID;
expoData.spiketimes.Times{1} = []; % Suppress this from analysis

% DEBUG PLOT - show the offset spiketimes 
if 0
  ok = ( spk_passID >= 1); %#ok<UNRCH>
  trigger_times = double(expoData.passes.StartTimes(1:2:end)) / 1e4;
  
  clf, hold on 
  plot( hekaData.Spikes.timeSeconds(spk_select & ok),   ...
        0*trigger_times(spk_passID(spk_select & ok)) +  ...
        time(order(hekaData.Spikes.timeIdx(spk_select & ok))), 'o')    
  plot( hekaData.Spikes.timeSeconds(~spk_select & ok),   ...
        0*trigger_times(spk_passID(~spk_select & ok)) +  ...
        time(order(hekaData.Spikes.timeIdx(~spk_select & ok))), 's')
end


