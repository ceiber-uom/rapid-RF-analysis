function [hekaData, expoData] = import_HEKA(varargin)
% import_HEKA: Load data and plot from HEKA .dat file with expo stimulus.
%   odat = load_expoHeka(filename, options) 
% 
% Update history
% 09-Jul-2016 CDE Wrote it, using code from /NANOEXPO/

%% Set up default options structure
options = utils.setup_options(varargin, ...
                             'extractSpikes',    1, ...
                             'removeLineNoise',  1, ... 
                             'saveDIR',          [pwd  filesep 'data' filesep], ...
                             'saveMAT',          1, ...
                             'saveXML',          1); 

%% Select recording if not specified
[heka_dat, expo_xml] = selectRecording(varargin);
if ~exist(heka_dat,'file'), return, end

disp(heka_dat)
all_hekaData = importHekaData(heka_dat);

name_of = @(n) regexp(n,'(?<=[/\\])[^/\\]+$','match','once');

if isempty(expo_xml)
    %%
    hekaData = all_hekaData;
    expoData = struct([]);    
    if isfield(options,'Protocol'), protocol = options.Protocol(1);
        if isnan(protocol) || isempty(protocol) || protocol > length(all_hekaData.Data)
            warning('importHekaData:badProtocolNumber' , ...
                    '%d is not a valid protocol number.', protocol)                
             return
        end
    else return % Protocol not specified
    end

    % Parse requested Protocol
    hekaData = selectHekaData(all_hekaData, protocol);

    if options.removeLineNoise
        hekaData.PassData = removeLineNoise_SpectrumEstimation(hekaData.PassData.', ...
                                          1/hekaData.SweepHeader.SampleInterval).';
    end

    if options.extractSpikes
        [hekaData, options] = detectSpikes(hekaData,options);
    end

    % Save if requested 
    if ~options.saveMAT, return, end

    mat_file = strrep(heka_dat,'.dat', sprintf('#%d[%s].mat', ...
                      protocol, hekaData.SweepHeader.EntryName));

    mat_file = regexprep(mat_file,'^.*[/\\](?=[^/\\]+$)',...
                         strrep(options.saveDIR,'\','\\'));
    if ~exist(options.saveDIR,'dir'), mkdir(options.saveDIR), end
    disp(['Saving ' mat_file]);
    save(mat_file,'expoData','hekaData','options');
    
else for ff = 1:length(expo_xml) % Load EXPO data
    %%
    if ~exist(expo_xml{ff},'file'), continue, end
    
    disp(expo_xml{ff})
    expoData = Tools.ex_ReadExpoXML_sc(expo_xml{ff}, 0);
    expoData.FileName = name_of(expo_xml{ff});
    hekaData = selectHekaData(all_hekaData, ff, expoData);
    
    if options.removeLineNoise
        hekaData.PassData = Tools.removeLineNoise_SpectrumEstimation(hekaData.PassData.', ...
                                          1/hekaData.SweepHeader.SampleInterval).';
    end

    if options.extractSpikes
        [hekaData, options] = detectSpikes(hekaData,options);
    end

    if options.extractSpikes && exist(expo_xml{ff},'file') && options.saveXML

        % import into expoData
        expoData.spiketimes.IDs = 0;
        expoData.spiketimes.Channels = 1;
        expoData.spiketimes.Times{1} = cast(hekaData.Spikes.timeSeconds .* ...
                                        expoData.environment.Conversion.TickDuration, ...
                                        'like', expoData.passes.StartTimes(1)) + ...
                                        expoData.passes.StartTimes(1);
        if ~exist(options.saveDIR,'dir'), mkdir(options.saveDIR), end
        addSpikesToXML(expo_xml{ff}, expoData, options.saveDIR);

    end

    if ~options.saveMAT && numel(expo_xml) == 1, return, end

    if exist(expo_xml{ff},'file'),
         mat_file = strrep(expo_xml{ff},'.xml','.mat');
    else mat_file = strrep(heka_dat,name_of(heka_dat), ... 
                          [expoData.FileName '.mat']);
    end

    mat_file = regexprep(mat_file,'^.*[/\\](?=[^/\\]+$)',...
                         strrep(options.saveDIR,'\','\\'));
    if ~exist(options.saveDIR,'dir'), mkdir(options.saveDIR), end
    disp(['Saving ' mat_file]);
    save(mat_file,'expoData','hekaData','options');
    
    end % <for>
end


if nargout == 0
    assignin('caller','expoData',expoData);
    assignin('caller','hekaData',hekaData);
    clear
end


end




function [f_dat, f_xml] = selectRecording(args)

persistent fpath
if isempty(fpath) || all(fpath == 0), fpath = '.\Heka'; end

args = args(cellfun(@ischar,args));
named = @(x) strncmpi(args,x,length(x));

find_xml = find(~cellfun(@isempty,regexp(args,'\.xml$')),1);
find_dat = find(~cellfun(@isempty,regexp(args,'\.dat$')),1);

if ~isempty(find_xml), f_xml = args{find_xml}; else f_xml = ''; end
if ~isempty(find_dat), f_dat = args{find_dat}; else f_dat = ''; end
    
if ~exist(f_dat,'file') && ~exist(f_xml,'file')
    [fname, fpath, fsel] = uigetfile({'*.dat', 'HEKA data file'; ...
                                      '*.xml', 'Expo XML file'}, [], ...
                                      fpath);
    if fsel == 0, return % No file selected
    elseif fsel == 1
        f_dat = [fpath fname];
    elseif fsel == 2     % .xml file selected  
        f_xml = [fpath fname];        
    end
    clear fname fsel
end

if ~exist(f_dat,'file')
    fp2 = regexp(f_xml,'^.*[/\\](?=[^/\\]+$)','match','once'); 
    f_dat = dir([fp2 '*.dat']);
    if length(f_dat) == 1
        f_dat = [fp2 f_dat.name];
    else
        [fname, fp2] = uigetfile({'*.dat', 'HEKA data file'}, ...
                           'Please select a HEKA data file' , fp2);
        f_dat = [fp2 fname];
    end
end


if ~exist(f_xml,'file') % Load (or auto-load) files
    f_xml = dir([fpath '*.xml']);

    if any(named('-all')) % Load all files
        fname = regexp(f_dat,'(?<=[/\\])[^/\\]+$','match','once');
        index = str2double(regexpi(fname,'(?<=Cell_?)\d+','match','once'));
        
        if isempty(index) || isnan(index)
            error('Could not determine Cell_[id] from %s', fname)
            % args(named('-all')) = [];
            % load_expoHeka(args{:});
            % f_dat = ''; return
        end
        
        RF_id = cellfun(@(n) max([NaN str2double(regexpi(n, ...
                           '(?<=RF_?)\d+','match','once'))]), {f_xml.name});
        if ~any(RF_id == index)
            error('No files for RF_%d found in the list %s', index, ...
                                           sprintf('\n%s',f_xml.name))
            % f_dat = ''; return
        end
        
        f_xml = {f_xml(RF_id == index).name};
        
        % Because "10" precedes "2" in filename order
        order = str2double(regexp(f_xml,'(?<=#)[\d]+','match','once'));
        [~,order] = sort(order,'ascend');
        
        f_xml = strcat(fpath,f_xml(order));
        
    elseif length(f_xml) == 1
        f_xml = [fpath f_xml.name];
    else
        [fname, fpath] = uigetfile({'*.xml', 'Expo XML file'}, ...
                            'Please select an Expo XML file' , fpath,'Multiselect','on');
        if ~ischar(fpath), f_xml = ''; return, end
        f_xml = strcat(fpath,fname);
    end
end

if ~iscell(f_xml), f_xml = {f_xml}; end

end

function [dat, options] = detectSpikes(dat, options)

if isfield(options,'spikes'), opts = options.spikes; 
else                          opts = struct;
end

if ~isfield(opts,'threshold'),   opts.threshold = -0.04;   end % -40mV threshold
if ~isfield(opts,'doWaveforms'), opts.doWaveforms = 1; end

[~,pidx] = findpeaks(dat.PassData(:), 'MinPeakHeight',opts.threshold);

[time,pass] = ind2sub(size(dat.PassData),pidx);

dat.Spikes.timeIdx = time';
dat.Spikes.passIdx = pass';

startTimes = [dat.PassInfo.Time];
spikeTimes = (time * dat.SweepHeader.SampleInterval) + startTimes(pass).';

dat.Spikes.timeSeconds = spikeTimes';

if opts.doWaveforms
    if ~isfield(opts,'waveSamples'), opts.waveSamples = -15:32; end
    idx = bsxfun(@plus,pidx, ones(size(pidx))*opts.waveSamples(:)');
    
    idx(idx < 1) = 1;
    idx(idx > numel(dat.PassData)) = numel(dat.PassData);
    
    dat.Spikes.spikeWaveforms = dat.PassData(idx).';
end

options.spikes = opts;

end

function dat = importHekaData(filename)

filename = strrep(filename,'.dat','.pgf');
fi = fopen(filename, 'rb','ieee-le');

% Base of STIM tree file
fi_name    = fread(fi, 4, 'uint8=>char'); %#ok<NASGU>
n_levels   = fread(fi, 1, '*int32');
level_size = fread(fi, double(n_levels), '*int32');

% Get the tree from the file
PGF_tree = readHekaTree(fi, level_size);
PGF_tree = translateHekaTree_PGF(PGF_tree);
fclose(fi);

filename = strrep(filename,'.pgf','.pul');
fi = fopen(filename, 'rb','ieee-le');

%% Base of PULSE tree file
fi_name    = fread(fi, 4, 'uint8=>char'); %#ok<NASGU>
n_levels   = fread(fi, 1, '*int32');
level_size = fread(fi, double(n_levels), '*int32');

% Get the tree from the file
PUL_tree = readHekaTree(fi, level_size);
PUL_tree = translateHekaTree_PUL(PUL_tree);
fclose(fi);


%%
filename = strrep(filename,'.pul','.dat');

fi = fopen(filename, 'rb','ieee-le');
wave = fread(fi,inf,'int16');
fclose(fi);

if isempty(wave), 
    error('importHekaData:problemReading: The data for %s could not be loaded into MATLAB.', filename) 
end

DAT_tree.wave = [];

for sweep = 1:length(PUL_tree.node)

    s = [PUL_tree.node(sweep).node(:).Data]/2;
    w = [PUL_tree.node(sweep).node(:).DataPoints];
    
    r = mean(diff(s))/mean(w);
    
    if r == round(r)
        for i = 1:length(s)
            DAT_tree(sweep).wave(:,i,:) = reshape(wave(s(i)+(1:w(i)*r)),w(i),r);     %#ok<AGROW>
        end
    else
        
        warning(['importHekaData:inconsistantData:' PGF_tree.node(sweep).EntryName],...
                'The data layout for sweep #%d[%s] could not be determined. wave returned as cell array.', ...
                sweep, PGF_tree.node(sweep).EntryName)

        if sweep == length(PUL_tree.node)
             s(end+1) = length(wave);
        else s(end+1) = PUL_tree.node(sweep+1).node(1).Data/2;
        end
        for i = 1:length(s)-1
            DAT_tree(sweep).wave{i} = wave(s(i)+(1:(s(i+1)-s(i))));
        end
    end
end

%% Format and return

PGF_tree.Protocol = PGF_tree.node;
PUL_tree.Protocol = PUL_tree.node;

dat.PGF = rmfield(PGF_tree,'node');
dat.PUL = rmfield(PUL_tree,'node');
dat.Data = DAT_tree;

end

function odat = selectHekaData(dat, id, expo)

if exist('expo','var') % Select using expo passData

    passes = cellfun(@numel,{dat.PUL.Protocol.node}); 
    passes = (passes == numel(expo.passes.IDs)/2);
    hashNo = str2double(regexp(expo.FileName,'(?<=#)[\d]+','match','once'));
    nFiles = evalin('caller','numel(expo_xml)');
    
    if length(unique({dat.PGF.Protocol(passes).EntryName})) > 1
                       
      % Spatial Frequency and Contrast both have 12 passes, resulting in 
      % bad confusion here!
        if ~isempty(strfind(lower(expo.FileName),'[sf')) % do SF!
            passes = passes & strcmp({dat.PGF.Protocol.EntryName},'SF');
        elseif ~isempty(strfind(lower(expo.FileName),'[contrast')) % do contrast!
            passes = passes & strcmp({dat.PGF.Protocol.EntryName},'Contrast');
        else
            arrayfun(@(d) fprintf('%3d | %s\n', d.NumberSweeps, d.EntryName), ...
                          dat.PGF.Protocol(passes))
            error('TODO')
        end
      % ON and OFF_size_ms use the same PGF Structure
        if ~isempty(strfind(lower(expo.FileName),'_size_ms')) % do SF!
            warning('importHekaData:size_multiselect' , ...
                   ['When importing "on_size_ms" and "off_size_ms", please select all' ... 
                    ' files to ensure correct import.'])
        end
    end
        
    if id == 1 && (sum(passes) ~= nFiles)
        files = evalin('caller','cellget(name_of,expo_xml)');
        warning('importHekaData:multipleMatches' , ...
                '%d files [%s%c] in HEKA match the signature(s) from %d XML files: %s\n', ...
                 sum(passes), sprintf('#%d,',find(passes)),8, numel(files),...
                 sprintf('\n%s',files{:}))
    end
    
    
    protocol = find(passes);
    valid_hash = ~(isempty(hashNo) || isnan(hashNo) || hashNo < 1 || ...
                           hashNo > length(dat.Data)); 

    
    
    if sum(passes) >= id, 
        protocol = protocol(id);
        fprintf('Protocol %d in HEKA matches the signature from\n%s.\n', ...
                 protocol, expo.FileName)
    elseif valid_hash && passes(hashNo)
        protocol = hashNo;
        fprintf('Protocol %d in HEKA matches the signature from\n%s.\n', ...
                 protocol, expo.FileName)
    elseif sum(passes) > 1
        protocol = protocol(1);
        warning('importHekaData:multipleMatches' , ...
                '[%s%c] in HEKA matches the signature from\n%s.\n', ...
                 sprintf('#%d,',find(passes)),8, expo.FileName)
    elseif ~valid_hash
        error('A valid protocol number could not be determined from %s', ...
                 expo.FileName)
    else
        error('[No files] in HEKA match the signature from\n%s.\n', ...
                 expo.FileName)
    end
    
    id = protocol; % Pick this one...
    clear passes valid_hash hashNo protocol files nFiles expo

end


odat.FileHeader    = rmfield(dat.PUL,'Protocol');
odat.SweepHeader   = dat.PGF.Protocol(id);
odat.RecordingInfo = rmfield(dat.PUL.Protocol(id),'node');
odat.PassInfo      = dat.PUL.Protocol(id).node;
odat.PassData      = dat.Data(id).wave;

V_unit = odat.PassInfo(1).DataVfactor; % Volts/ADC-unit
A_unit = odat.PassInfo(1).DataAfactor; % Amps/ADC-unit
N_dims = size(odat.PassData,3);

if iscell(odat.PassData) % massive assumption here!!!    
    odat.PassData = cellfun(@(s)V_unit*s, odat.PassData,'UniformOutput',false);    
    return
end

% Convert to real units
odat.PassData(:,:,1) = odat.PassData(:,:,1)*V_unit;
if odat.PassInfo(1).Leak && N_dims > 1 % second trace = leak
   odat.PassData(:,:,2) = odat.PassData(:,:,2)*V_unit;
end
if odat.PassInfo(1).SecondTrace && N_dims > 1 % last trace = current
   odat.PassData(:,:,end) = odat.PassData(:,:,end)*A_unit;
end

for ii = length(odat.PassInfo):-1:1, 
    odat.PassInfo(ii).Time = odat.PassInfo(ii).Time - odat.PassInfo(1).Time;
end





end

function addSpikesToXML(name, dat, folder)

% read XML as string
fi = fopen(name,'r');
src = fread(fi,inf,'*char')';
fclose(fi);

% flip <Events .. Spike="0" to ="1">
idx = regexp(src,'(?<=Events [^/]*)Spike="'); 

for ii = 1:length(idx)
    % strangely, couldn't do this in regex string??
    if isempty(regexp(src(idx(ii)-1),'\s', 'once')), continue, end
    src(idx(ii)+7) = '1';
end

% add actual spiketime data
idx = regexp(src,'</ExpoXData>');
spk_str = regexprep(num2str(dat.spiketimes.Times{1}),'\s+',',');
spk_str = sprintf(['\t<Spikes Total="%d">\n\t' ... 
                   '\t<Spike ID="0" Channel="1" Times="%s"/>\n' ...
                   '\t</Spikes>\n'], length(dat.spiketimes.Times{1}), spk_str);
src = [src(1:idx-1) spk_str src(idx:end)];

% write to new file
name = strrep(name,'.xml','_spike.xml');
name = regexprep(name,'^.*[/\\](?=[^/\\]+$)',strrep(folder,'\','\\'));
fi = fopen(name,'w');     disp(['Saving ' name]);
fwrite(fi,src);
fclose(fi);

               
end







%% Low-level functions: Recursive tree read function
function tree = readHekaTree(fi, levels)

tree.raw = fread(fi,levels(1),'*uint8');
n_children = fread(fi,1,'*int32');

if length(levels) == 1 || n_children == 0, return, end

tree.node = readHekaTree(fi, levels(2:end));
for ii = 2:n_children
    tree.node(ii) = readHekaTree(fi, levels(2:end));
end

end

% Recursive tree translate functions.
function tree = translateHekaTree_PGF(tree)

data_length = length(tree.raw);

switch data_length
    case 2   %   V8.8 PGF Header        
        tree.Version = typecast(tree.raw,'int16');
    
    case 408 %   V8.8 PGF Stimulus Record. One per Expo XML
        tree.FileName         =     char(tree.raw(  1:14 ))';         % Source File
        tree.EntryName        =     char(tree.raw( 15:28 ))';         % Identifier
        tree.SampleInterval   = typecast(tree.raw( 29:36 ),'double'); % Seconds
        tree.FilterFactor     = typecast(tree.raw( 37:44 ),'double'); % oversampling factor
        tree.SweepInterval    = typecast(tree.raw( 45:52 ),'double'); % Repetition-Interval
        tree.NumberSweeps     = typecast(tree.raw( 53:56 ),'int32');  % Number of sweeps in Series
        tree.NumberRepeats    = typecast(tree.raw( 57:60 ),'int32');  % Number of Seq. repeats
        tree.RepeatWait       = typecast(tree.raw( 61:68 ),'double'); % Wait between repeats
        tree.LinkedSequence   =     char(tree.raw( 69:82 ))';
        tree.LinkedWait       = typecast(tree.raw( 83:90 ),'double'); % Wait between linked sequences
        tree.LeakCount        = typecast(tree.raw( 91:94 ),'int32');  % number of leak sweeps
        tree.LeakSize         = typecast(tree.raw( 95:102),'double'); % rel. amplitude of leak pulse
        tree.LeakHolding      = typecast(tree.raw(103:110),'double');
        tree.LeakAlternate    =  logical(tree.raw(  111  ));          % normal or alt leak protocol
        tree.AltLeakAveraging =  logical(tree.raw(  112  ));          % normal or alt leak while averaging
        tree.LeakDelay        = typecast(tree.raw(113:120),'double'); % Seconds

        for t = 1:3
            tree.Trig(t).TriggerSegment   = typecast(tree.raw(t*28+92+( 1:2)),'int16');
            tree.Trig(t).TriggerTime      = typecast(tree.raw(t*28+92+( 3:10)),'double'); % Seconds, within trig Segment
            tree.Trig(t).TriggerLength    = typecast(tree.raw(t*28+92+(11:18)),'double'); % Length of Trig Pulse, sec
            tree.Trig(t).TriggerAmplitude = typecast(tree.raw(t*28+92+(19:26)),'double'); % Amplitude of Trig Pulse 
            tree.Trig(t).TriggerDac       = typecast(tree.raw(t*28+92+(27:28)),'int16');
        end
        
        tree.NumberOfTriggers = typecast(tree.raw(205:206),'int16');  % Number actually used; 0<x<4
        tree.RelevantXSegment = typecast(tree.raw(207:208),'int16');  % usually that which changes
        tree.RelevantYSegment = typecast(tree.raw(209:210),'int16');  % usually that which changes        
        
        enum = {'0:WriteEnabled','1:WriteDisabled','2:NoWriteNoShow','3:WriteButNoShow','+:Unkown'};        
        tree.WriteMode        = enum{min(tree.raw(211)+1, length(enum))};  
        enum = {'0:ModeInc','1:ModeDec','2:ModeIncInterleaved','3:ModeDecInterleaved','4:ModeAlternate','+:Unkown'};
        tree.IncrementMode    = enum{min(tree.raw(212)+1, length(enum))};                  
        tree.TotalSweepLength = typecast(tree.raw(213:216),'int32');  % Total sweeplength, in units of sample intervals
        tree.MaxSweepLength   = typecast(tree.raw(217:220),'int32');  % Max length of sweep to be shown.
        tree.InputChannels    = typecast(tree.raw(221:222),'int16');  % Number of input channels. Default 1
        
        enum = {'0:NoGUpdate','1:SwGSlow','2:SwGFast','3:SwGBoth','4:SeGSlow','5:SeGFast','6:SeGBoth','+:Unkown'};
        tree.GUpdate        = enum{min(tree.raw( 223 )+1, length(enum))};
        tree.RelAbsPot      =  logical(tree.raw( 224 ));            % absolute or relative potentials
        tree.HasContinuous  =  logical(tree.raw( 225 ));
        tree.LogIncrement   =  logical(tree.raw( 226 ));        
        tree.StimDac        = typecast(tree.raw(227:228),'int16');
        tree.Adc_1          = typecast(tree.raw(229:230),'int16');
        tree.Adc_2          = typecast(tree.raw(231:232),'int16');
        tree.YUnit_1        =     char(tree.raw(233:234))';
        tree.YUnit_2        =     char(tree.raw(235:236))';  
        tree.VmembIncrement = typecast(tree.raw(237:240),'single');  % Number of Seq. repeats

        enum = {'0:TrigNone','1:TrigSeries','2:TrigSweep','+:Unkown'};
        tree.ExtTrigger     = enum{min(tree.raw( 241 )+1, length(enum))};
        tree.FileTemplate   =  logical(tree.raw( 242 ));
        
        tree.StimKind.LockinActive = logical(bitand(tree.raw( 243 ),1,'uint8'));
        tree.StimKind.FuraActive   = logical(bitand(tree.raw( 243 ),2,'uint8'));
        tree.StimKind.useScanRate  = logical(       tree.raw( 244 )           ); % technically bit 15 ... 
        tree.LockInCycle     = typecast(tree.raw(245:252),'double');
        tree.LockInAmplitude = typecast(tree.raw(253:260),'double');
        
        
        enum = {'0:VmembSweepIncr','1:VmembValue','2:VmembSeriesIncr','+:Unkown'};
        tree.FuraOn           =  logical(tree.raw( 261 )); 
        tree.VmembMode        = enum{min(tree.raw( 262 )+1, length(enum))};
        tree.FuraTotLength    = typecast(tree.raw(263:270),'double');
        tree.FuraDelay        = typecast(tree.raw(271:278),'double');
        tree.FuraLength_1     = typecast(tree.raw(279:286),'double');
        tree.FuraLength_2     = typecast(tree.raw(287:294),'double');
        tree.FuraWaveLength_0 = typecast(tree.raw(295:302),'double');
        tree.FuraWaveLength_1 = typecast(tree.raw(303:310),'double');
        tree.FuraWaveLength_2 = typecast(tree.raw(311:318),'double');
        tree.FuraRepeats      = typecast(tree.raw(319:320),'int16');
        tree.LockInSkip       = typecast(tree.raw(321:324),'int32');
        tree.LockInVReversal  = typecast(tree.raw(325:332),'double');
        
        enum = {'0:Loff','1:Lnormal','2:Lpwlinear','+:Unkown'};        
        tree.LockInMode  = enum{min(tree.raw( 333 )+1, length(enum))};
        tree.LockInShow  =  logical(tree.raw( 334 ));
        tree.ConfigMacro =     char(tree.raw(335:350))';
        tree.EndMacro    =     char(tree.raw(351:366))';

        enum = {'0:AllAmplModes','1:VCAmplMode','2:CCAmplMode','+:Unkown'};        
        tree.AmplModeKind      = enum{min(tree.raw( 367 )+1, length(enum))};
        tree.NoStartWait       =  logical(tree.raw( 368 ));
        tree.ActualInChannels  = typecast(tree.raw(369:370),'int16');
        tree.ActualOutChannels = typecast(tree.raw(371:372),'int16');
        tree.ActualAdc_1       = typecast(tree.raw(373:374),'int16');
        tree.ActualAdc_2       = typecast(tree.raw(375:376),'int16');
        % tree.RawUnprocessed  = tree.raw(377:end);
    
    case 50 %   V8.8 PGF Stimulus Segment.

        enum = {'0:SegmentConstant',  '1:SegmentRamp',     '2:SegmentConditioning',...
                '3:SegmentContinuous','4:SegmentConstSine','5:SegmentRampSine','+:Unkown'};
        tree.Class           =     enum{tree.raw(  1  )+1}; % constant, ramp, or condit
        tree.IsHolding       =  logical(tree.raw(  2  ));   % TRUE if Voltage is to be set
        tree.Voltage         = typecast(tree.raw( 3:10),'double');  % to holding
        tree.Duration        = typecast(tree.raw(11:18),'double');
        tree.DeltaVFactor    = typecast(tree.raw(19:26),'double');
        tree.DeltaVIncrement = typecast(tree.raw(27:34),'double');
        tree.DeltaTFactor    = typecast(tree.raw(35:42),'double');
        tree.DeltaTIncrement = typecast(tree.raw(43:50),'double');

    otherwise
        warning('importHekaData:unknownPGLformat','PGF record of unknown length: %d bytes', data_length)
        tree.Raw = char(tree.raw)';
end

% Clean up text fields
fn = fieldnames(tree);
for ii = 1:length(fn), if ischar(tree.(fn{ii})),  %#ok<ALIGN>
    tree.(fn{ii}) = tree.(fn{ii})(1:min([find(tree.(fn{ii})==0,1)-1 end]));
end, end

tree = rmfield(tree,'raw');
if ~isfield(tree,'node'), return, end

node = translateHekaTree_PGF(tree.node(1));
for ii = 2:length(tree.node), node(ii) = translateHekaTree_PGF(tree.node(ii)); end
tree.node = node;


end
function tree = translateHekaTree_PUL(tree)

data_length = length(tree.raw);

switch data_length
        
        
    case 438 % PUL Header, but I'm not confident that this is 100% stable.
        
        tree.Version = typecast(tree.raw(1:2),'int16');        
        tree.Comments = char(tree.raw(3:end-8))';
        tree.StartTime = typecast(tree.raw(end-7:end),'double'); % oversampling factor
        
        % Also grab the singleton node beneath!
        tree.raw  = tree.node.raw;
        tree.node = tree.node.node;
        % A GroupRecord describes a group of series, such as patch and
        % whole cell currents obtained simultanuously, or groups of series
        % obtained in sequence under different sets of conditions
        tree.RecordLabel = char(tree.raw( 1:14))';
        tree.Text        = char(tree.raw(15:94))';
        tree.ExperimentNumber = typecast(tree.raw(95:98),'int32');
        tree.ExtraLongReal    = typecast(tree.raw(99:106),'double');   
                
    case 440

        % Patch parameters, usually obtained before series        
        tree.RecordTime        = typecast(tree.raw(  1:8  ),'double'); % Seconds
        tree.Bandwidth         = typecast(tree.raw(  9:16 ),'double'); % Hertz
        tree.PipettePotential  = typecast(tree.raw( 17:24 ),'double'); % Volts
        tree.CellPotential     = typecast(tree.raw( 25:32 ),'double'); % Volts
        tree.PipetteResistance = typecast(tree.raw( 33:40 ),'double'); % Ohms
        tree.SealResistance    = typecast(tree.raw( 41:48 ),'double'); % Ohms
        tree.BackgroundNoise   = typecast(tree.raw( 49:56 ),'double'); % Amperes, RMS
        tree.Temperature       = typecast(tree.raw( 57:64 ),'double'); % Degrees C
        tree.PipettePressure   = typecast(tree.raw( 65:72 ),'double'); % in cm H20
        
        tree.UserParam.Name    =     char(tree.raw( 81:94 ))';
        tree.UserParam.Value   = typecast(tree.raw( 73:80 ),'double');
        tree.UserParam.Unit    =     char(tree.raw( 95:96 ))';
        tree.UserParam(2).Value = typecast(tree.raw( 97:104 ),'double');
        tree.UserParam(2).Name  =     char(tree.raw(105:118))';
        tree.UserParam(2).Unit  =     char(tree.raw(119:120))';
        
        enum = {'0:InOut','1:OnCell','2:OutOut','3:WholeCell','4:PCClamp', ...
                '5:VClamp','6:NoRec','7:TestInt','8:TestExt','+:Unkown'};
        tree.RecordingMode = enum{min(tree.raw(121)+1, length(enum))};
        tree.VideoMode = tree.raw( 122 );
        tree.Comment = char(tree.raw(123:202))';

        % EPC9State : EPC9StateType;
        % InternalSolution,ExternalSolution : INT32; (* solutions according to solutioncode as stored in '*.sol' *)
        % ExtraYUnit1 : ARRAY [0..1] OF CHAR;
        % ExtraYUnit2 : ARRAY [0..1] OF CHAR;
        % DispYUnit1 : ARRAY [0..3] OF CHAR; (* reserved by PULSE *)
        % DispYUnit2 : ARRAY [0..3] OF CHAR; (* reserved by PULSE *)
        % FuraK : LONGREAL;
        % FuraMin : LONGREAL;
        % FuraMax : LONGREAL;
        % LockInExtPhase : LONGREAL;
        % Timer : LONGREAL; (* Reserved for extensions *) 
        % ExtCalPhase : LONGREAL;
        % ExtCalAttenuation : LONGREAL;
        % PLPhase : LONGREAL;
        % PLPhaseY1 : LONGREAL;
        % PLPhaseY2 : LONGREAL;
        % ExtCalValid : BOOLEAN;
        % PLPhaseValid : BOOLEAN;
        
        % (* EPC9 state record *)EPC9State : EPC9StateType ????? 
        % EPC9StateType = ARRAY[0..103] OF BYTE; (* E9Panel.StateType *)
    
    case 216
        
        tree.Time          = typecast(tree.raw( 1:8 ),'double'); % Seconds (unix)        
        tree.StimCount     = typecast(tree.raw( 9:12),'int32'); % relevant entry in StimTree
        tree.SweepCount    = typecast(tree.raw(13:16),'int32'); % .. at time of acquisition        
        tree.AverageCount  = typecast(tree.raw(17:20),'int32'); % number of on-line averages
        tree.Leak          =  logical(tree.raw( 21 )); % TRUE, if Record has a leaksweep
        tree.SecondTrace   =  logical(tree.raw( 22 )); % TRUE, if Record has a 2nd trace        
        tree.Label         =     char(tree.raw(23:36))';
        tree.DataPoints    = typecast(tree.raw( 37:40 ),'int32'); % Number of data points in RAM
        tree.Data          = typecast(tree.raw( 41:44 ),'int32'); % Offset in raw data file
        tree.DataPointer   = typecast(tree.raw( 45:48 ),'uint32'); % ... Raw data ADDRESS in memory, 
        
        % EPC 7/9 settings
        tree.DataAfactor   = typecast(tree.raw( 49:56 ),'double'); % Amperes/ADC-unit
        tree.DataVfactor   = typecast(tree.raw( 57:64 ),'double'); %  Volts/ADC-unit
        tree.CSlow         = typecast(tree.raw( 65:72 ),'double'); % Capacitance, Farad
        tree.GSeries       = typecast(tree.raw( 73:80 ),'double'); % Series conductance, Siemens
        tree.RsValue       = typecast(tree.raw( 81:88 ),'double'); % Series resistance setting, Ohms
        tree.Mconductance  = typecast(tree.raw( 89:96 ),'double'); % Results of pre-analysis
        tree.ZeroCurrent   = typecast(tree.raw( 97:104),'double'); % Amperes
        tree.OnlineXResult = typecast(tree.raw(113:120),'double'); % Param determined by last analysis
        tree.OnlineYResult = typecast(tree.raw(105:112),'double'); % Param determined by last analysis        
        tree.TotalPoints   = typecast(tree.raw(121:124),'int32'); % Total data points in file
        tree.Offset        = typecast(tree.raw(125:128),'int32'); % Offset of loaded data
        tree.SweepKind     = typecast(tree.raw(129:130),'int16'); % (* meaning of bits:ExtendedBit => extended fields definedFuraBit => Fura.ActiveLockInBit => LockIn.ActiveDataFormatBit => "DataFormat" field validLittleEndianBit => byte sequence"Mac" = cleared"DOS" = setScanRateBit => set: use scan-rate!*)
        tree.FuraPoints    = typecast(tree.raw(131:134),'int32'); % Number of data points in RAM
        tree.FuraData      = typecast(tree.raw(135:138),'int32'); % Offset in raw data file
        tree.FuraPointer   = typecast(tree.raw(139:142),'uint32'); % Raw data ADDRESS in memory        
        tree.OnlineYResult(2) = typecast(tree.raw(143:150),'double'); % Param determined by last analysis
        tree.OnlineXResult(2) = typecast(tree.raw(151:158),'double'); % Param determined by last analysis
        tree.DispFactor(1) = typecast(tree.raw(159:166),'double'); % reserved by PULSE
        tree.DispFactor(2) = typecast(tree.raw(167:174),'double'); % reserved by PULSE
        
        enum = {'int16','int32','single','double'};
        tree.DataFormat    = enum{tree.raw(175)+1};
        
        enum = {'0:Time','1:Amplitude','2:LnAmplitude','3:LogAmplitude', ...
                         '4:Frequency','5:LnFrequency','6:LogFrequency','+:Unkown'};
        tree.DataAbscissa  = enum{min(tree.raw(176)+1,length(enum))};        
        tree.Timer         = typecast(tree.raw(177:184),'double'); % Param determined by last analysis
        tree.FuraData      = typecast(tree.raw(185:188),'int32'); % Offset in raw data file
        tree.VideoTime     = typecast(tree.raw(189:196),'double'); % Param determined by last analysis

        tree.Spares = char(tree.raw(197:end))';
        
    otherwise
        warning('importHekaData:unknownFormat','Record of unknown length: %d bytes', data_length)
        tree.Raw = char(tree.raw)';
end

% Clean up text fields
fn = fieldnames(tree);
for ii = 1:length(fn), if ischar(tree.(fn{ii})),  %#ok<ALIGN>
    tree.(fn{ii}) = tree.(fn{ii})(1:min([find(tree.(fn{ii})==0,1)-1 end]));
end, end
    
tree = rmfield(tree,'raw');
if ~isfield(tree,'node'), return, end

node = translateHekaTree_PUL(tree.node(1));
for ii = 2:length(tree.node), node(ii) = translateHekaTree_PUL(tree.node(ii)); end
tree.node = node;


end

