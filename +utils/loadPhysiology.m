
function loadPhysiology(varargin)

named = @(n) strncmpi(varargin,n,length(n));

persistent fp
if isempty(fp) || any(fp == 0), 
    fp = regexprep(fileparts(mfilename('fullpath')),'([/\\])\+.*','$1MAT$1');
end

fn = cellfun(@(v) ischar(v) && any(v=='#'), varargin); 

if any(fn), fn = varargin{find(fn,1)}; 
    
    if exist(fn,'file'), disp(['Loading ' fn]), load(fn);
    elseif exist([fp fn],'file'), disp(['Loading ' fn]), load([fp fn]); 
    else
        % Parse code of the form "8.6.1 #9" 
        % also can accept "20180806_Cell_1 #9"
        dc_code = str2double(regexp(fn,'\d+','match'));
        
        if numel(dc_code) == 3 && dc_code(1) > 2e6                
            dc_code = [mod(dc_code(1),[1e4 1e2]) dc_code(2:end)];        
            dc_code(1) = (dc_code(1) - dc_code(2))/100; 
        elseif numel(dc_code) ~= 4, error('%s not a valid code string.', fn), 
        end
        fn = sprintf('%02d%02d_Cell_%d',dc_code(1:3));
        
        list = dir([fp '*'  fn '*.mat']);
        if isempty(list), error('%s not found in %s', fn, fp), end
        hash_no = cellfun(@(s) str2double(regexp(s,'(?<=#\w*)\d+','match')), {list.name}); 
        
        if ~any(hash_no == dc_code(end)), 
            error('%s #%d not found in [%s%s].', fn, dc_code(end), sprintf('#%d,',hash_no),8)
        end

        fn = list(hash_no == dc_code(end)).name;
        disp(['Loading ' fn]), load([fp fn]);
    end
else
    [fn,fp] = uigetfile('*.mat','Select data',fp);    
    disp(['Loading ' fn])
    load([fp fn]); 
end

%%
fs = 1./hekaData.SweepHeader.SampleInterval;            %#ok<NODEF>
passes = double(expoData.passes.BlockIDs(2:2:end));
time = (0:size(hekaData.PassData,1)-1) / fs ;

delay = diff(expoData.passes.StartTimes);
delay = mean(delay(1:2:end)) / 1e4; 
time = time - delay; 

duration = double(expoData.passes.EndTimes - ...
                  expoData.passes.StartTimes) / 1e4;
duration = mean(duration(2:2:end));

nPasses = length(passes);

%% Remove Action Potentials from membrane current
if ~any(named('-spike')), 
    hekaData = Tools.removeActionPotentials(hekaData); 
end

%% If needed, apply a time offset (if response has crept into next pass)
if isfield(options,'apply_offset')
    
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
end


%% Direct plot of membrane potential    

stim_rate = expoData.passes.events{2}.Data{1}(3);
stim_bar = @(i,h) [(i-1)/stim_rate ylim*[h;1-h] 1/stim_rate/2 ylim*[-1;1]*abs(h/3)];
nStimuli = round(duration * stim_rate);  

if any(named('-plot')), 
    
    figure(1), clf
    set(gcf,'Color','w','Position',[50 560 960 330],'Name','Membrane Potential')

    plot(time, 1000*hekaData.rawPassData, 'Color', [0 0 0 0.05])
    xlim([min(time) max(time)])

    for ss = 1:nStimuli,
        rectangle('Position',stim_bar(ss,0.1),'FaceColor',[0 0 0 0.5], 'EdgeColor','none')
        % rectangle('Position',stim_bar(ss+0.5),'FaceColor',[0  0  0  .4], 'EdgeColor','none')
    end

    set(gca,'XTickMode','manual','XTickLabel',strcat(get(gca,'XTickLabel'),' s'))
    set(gca,'YTickMode','manual','YTickLabel',strcat(get(gca,'YTickLabel'),' mV'))
    set(gca,'TickDir','out','box','off','Color','none')
end


%% Dump the loaded data into the caller workspace

assignin('caller','filename',strrep(fn,'_','\_'))
assignin('caller','hekaData',hekaData)
assignin('caller','expoData',expoData)

assignin('caller','time',time)
assignin('caller','passes',passes)
assignin('caller','nPasses',nPasses)
assignin('caller','nStimuli',nStimuli)
assignin('caller','stim_bar',stim_bar)

%% Get the anatomy if possible
if any(named('-anatomy')), 
    
    f = figure; pause(0.05)
    Tools.loadAnatomy('-check');
    im = get(f.Children,'Children');    
    
    if isempty(im)
        anatomy.plotCell = @(h)[]; 
    else
        
        tracing = fullfile('Anatomy','Filled Anatomy',strrep(im.UserData,'.png','.mat'));
        if exist(tracing,'file'), anatomy = load(tracing);
        else                      anatomy = struct;
        end
            
        anatomy.img = im.CData; 
        anatomy.X = im.XData;
        anatomy.Y = im.YData;
        anatomy.alpha = im.AlphaData;               
        anatomy.plotCell = @(h) imagesc(anatomy.X, anatomy.Y, anatomy.img,'AlphaData', anatomy.alpha,'Parent',h);
    end
    assignin('caller','anatomy',anatomy)

    close(f)
end

%% Compute PSTH (if requested)


if any(named('-psth'))

    wave_time = time; % Possible "jumps" - ignore for analysis 
    
    if any(named('-bin'))        
         time = varargin{find(named('-bin'))+1}; % e.g. -nK 2
    elseif isfield(options,'apply_bin_size')
         time = options.apply_bin_size; 
         fprintf('Using non-standard bin size %d (default: 100) for this file\n', time)        
    else time = 100;
    end
    
    time = 0:time:size(hekaData.PassData,1);
    time = time(1:end-1)/2 + time(2:end)/2;

    wave = zeros(numel(time), size(hekaData.PassData,2));

    for ii = 1:size(hekaData.PassData,1) % make PSTH 
        pass = hekaData.Spikes.passIdx == ii;
        if ~any(pass), continue
        else wave(:,ii) = hist(hekaData.Spikes.timeIdx(pass),time);
        end
    end

    time = wave_time(time); % If there's a jump, need to know
    % time = (time * hekaData.SweepHeader.SampleInterval) - delay;    

    psth.time = time;
    psth.wave = wave / median(diff(time));

    assignin('caller','psth',psth)

end

%% PCA or NNMF analysis (if requested)

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

    % Debug: spectrum of the stimulus S
    % clf
    % hz = linspace(0,1,sum(roi))*fs; 
    %  
    % loglog(hz(ok),abs(S(ok)))
    % axis tight, hold on
    % if ~exist('pp','var'), [~,pp] = max(abs(f1)); end
    % loglog(hz(ok),abs(F(ok,pp)))
    % loglog(hz(ok),abs(F(ok,pp) .* S(ok)'))
    % plot(hz(ok),abs(f1(pp))*ones(size(ok)),'--')

    % % Debug: plot max of ifft wave vs computed f1
    % clf
    % plot(max(f1_wave),abs(f1),'.')
    % hold on, plot(xlim,xlim,'Color',[0 0 0 0.3])

    di = mean(exp(1i*angle(f1))); 
    f1 = f1 / (di/abs(di)) * -sign(real(di)); 
    f0 = mean(hekaData.PassData(roi,:));

    assignin('caller','f0',f0)
    assignin('caller','f1',f1)
    assignin('caller','f1_wave',f1_wave);
end


if any(named('-nK')), 
    if any(named('-nK:')) % As string argument e.g. -nK:2
         nK = str2double(regexp(varargin{find(named('-nK'),1)},'\d+','match','once'));
    else nK = varargin{find(named('-nK'))+1}; % e.g. -nK 2
    end
else     nK = 12; 
end

if any(named('-nnmf'))
    
    if any(named('-psth')), 
         [score,coeff] = nnmf(max(psth.wave,0),nK);
    else [score,coeff] = nnmf(-min(1e3*hekaData.PassData,0),nK);
          score = -score;
    end
    
    coeff = coeff';  
    
    resting = 0; 
    
elseif any(named('-pca'))
    
    if any(named('-psth'))
        if ~any(time < 0), resting = 0; 
        else resting = median(reshape(psth.wave(time < 0,:),[],1));
        end
        [coeff,score] = pca((psth.wave - resting),'Centered',false);
    else
        if ~any(time < 0), resting =quantile(hekaData.PassData(:),0.2); 
        else resting = median(reshape(hekaData.PassData(time < 0,:),[],1));
        end        
        [coeff,score] = pca(1e3*(hekaData.PassData - resting),'Centered',false);
    end
    

else return
end

coeff = coeff(:,1:nK);  coe_scale = max(coeff);
score = score(:,1:nK);

coeff = coeff * diag(1./coe_scale);
score = score * diag(coe_scale);

assignin('caller','coeff',coeff)
assignin('caller','score',score)
assignin('caller','resting',resting)
assignin('caller','coe_scale',coe_scale)
    
end
    
    
    

