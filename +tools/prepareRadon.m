function dat = prepareRadon(dat, varargin)
% rdat = tools.prepareRadon(data, ... )
% rdat = tools.prepareRadon(data, activation_maps, ... )
% rdat = tools.prepareRadon(activations, [waves, time], ... )
% 
% Convert the input data into the correct format for analysis.inverseRadon
% 
% The output has the following fields: 
%  .ori [nBlocks x 1] - angle of the bar in degrees
%  .x   [nBlocks x 1] - distance of the bar to the centre in µm or deg
%  .y_all [nBlocks x nY] - matrix of responses to stimulation 
%  .time [1 x nT]  (if do_profiles) - time vector for waves
%  .wave [nT x nY] (if do_profiles) - component waveforms
% 
% Options for conversion:
% -no-wave : do not return .time and .wave, even if you find them. 
% -pixel   : set pixel conversion size, default 2.342818332 (Protti lab)
%            -LGN (Martin lab) is equivalent to -pixel 1
% -expo [expoData] : use supplied stimulus info (can auto-detect).
% -time [time]     : use supplied time vector (can auto-detect).
% -smooth          : do not apply "histo stairs vector" effect to
%                     psth-based waves. 
% -append : add the fields to 'dat' rather than returningd a new structure.
%           Example: dat = tools.prepareRadon( dat, '-append' )
% 
% if DATA is a struct, one of the following fields is selected 
%    (in descending order of preference): 
% 
%   .activations (e.g. from tools.load( ..., 'pca' or 'nnmf')
%   .f0, .f1     (e.g. from tools.load( ..., '-f01')
%   .psth.wave   (e.g. from tools.load( ..., '-psth')
%   .hekaData.PassData (if not found, not recommended)
% 
% if DATA is a struct and activation_maps is supplied, that is used instead
% if DATA.response_waves is found, .time and .wave are returned as well. 
% 
% v2.0 - 28 August 2022 - Calvin Eiber <ceiber@ieee.org>
%        Extracted from generic conversion code


named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

do_profiles = true;
if any(named('-no-wave')), do_profiles = false; end
if any(named('-pixel')),   um_per_pixel = get_('-pixel'); 
elseif any(named('-lgn')), um_per_pixel = 1; % degrees of visual field
else                       um_per_pixel = 2.342818332; % Screen conversion
end

make_time_histoStairsVector = false; 

if isstruct(dat)
  
  if all(isfield(dat,{'x','ori','y_all'})), return, end % nothing to change
 
  %% This is probably RAW data as loaded by tools.load

  if any(named('-expo')), expoData = get_('-expo');
  elseif isfield(dat,'expoData'), expoData = dat.expoData; 
  elseif nargin > 1 && isstruct(varargin{1}) && isfield(varargin{1},'ExpoVersion')
                              expoData = varargin{1};
  else error('please supply -expo expoData');
  end
  
  if nargin > 1 && isnumeric(varargin{1}) % dat, map syntax
                                     activations = varargin{1}; 
  elseif isfield(dat,'activations'), activations = dat.activations; 
  elseif isfield(dat,'f0'),   activations = [dat.f0; real(dat.f1); ...
                                                     imag(dat.f1)]';
  elseif isfield(dat,'psth'), activations = dat.psth.wave'; 
      % warning('radon:prepare:rawSpike','preparing to conduct radon analysis of PSTH')
  elseif isfield(dat,'hekaData'), activations = dat.hekaData.PassData';
      % warning('radon:prepare:rawWave','preparing to conduct radon analysis of raw Vm')
  else error('unable to locate input data.')
  end

  if numel(expoData) > 1 % tools.load (multiple files together)

      the_radon = contains(lower({expoData.FileName}),'radon');
      if any(named('-prot')), the_radon = get_('-prot');
      elseif sum(the_radon) == 1, the_radon = find(the_radon);
      else
          for ii = 1:numel(expoData)
              fprintf('%2d: %s\n', ii, expoData(ii).FileName)
          end
          error('multiple radon expoData: please set -protocol ID ')
      end

      activations = activations(dat.protocol == the_radon,:);
      expoData = expoData(the_radon);
  end

  if do_profiles
    if isfield(dat,'response_waves'), wave = dat.response_waves;
       if isfield(dat,'psth'), time = dat.psth.time;
                               make_time_histoStairsVector = true; 
       else                    time = dat.time;
       end
    else do_profiles = false; 
    end
  end

  
else
  %% Parse (wave,profile) input arguments

  if nargin > 1 && isnumeric(varargin{1}) % [waves, [time]]
    
    
    wave = varargin{1};
    
    infer_time = false; 

    if any(named('-time')), time = get_('-time'); 
    elseif nargin > 2 && isnumeric(varargin{2}), time = varargin{2};
    elseif evalin('caller','exist("time","var")')
        time = evalin('caller','time'); infer_time = true; 
    else do_profiles = false; time = []; 
    end

    if do_profiles && (length(time) ~= size(wave,1) || any(named('-psth')))
        if infer_time, time = evalin('caller','psth.time'); end     
        make_time_histoStairsVector = true; 
    end
  end
    
  activations = dat;
  if any(named('-expo')), expoData = get_('-expo');
  elseif evalin('caller','exist("expoData","expoData")')
      expoData = evalin('caller','expoData');
  else error('please supply -expo expoData');   
  end
  
end

%% From the selected data, construct the OUT vector

out = struct;   

% Extract Stimulus information
A = cellfun(@(p) p.Data{3}(1), expoData.passes.events(2:2:end));
X = cellfun(@(p) p.Data{3}(2), expoData.passes.events(2:2:end));
X(A >= 180) = -X(A >= 180); 
A(A >= 180) =  A(A >= 180) - 180; 

[out.ori,~,idx] = unique([A' X'],'rows');

out.x   = out.ori(:,2) * um_per_pixel; % um per pixel
out.ori = out.ori(:,1);

for ii = 1:length(out.ori)    
  out.y_all(ii,:) = mean(activations(idx == ii,:), 1);
end

if make_time_histoStairsVector && ~any(named('-smooth'))
    % convert to PSTH bins
    time = [reshape(time,1,[]) - mean(diff(time))/2 ...
            time(end) + mean(diff(time))/2];
    time = reshape([1;1] * time,1,[]);
    wave = cat(1, permute(wave,[3 1 2]), permute(wave,[3 1 2]));
    wave = reshape(wave,[],size(wave,3));
    wave = [0*wave(1,:); wave; 0*wave(1,:)];
end

if do_profiles
  out.time = time; 
  out.wave = wave;
end
clear expoData A X activation

if any(named('-append'))
  for f = fieldnames(out)', dat.(f{1}) = out.(f{1}); end
else dat = out; 
end

return
