
function result = linear(data, varargin)
%  Y = analysis.linear( X, [nK], [mode], ... ) 
% 
% Linear decomposition of the input matrix using PCA, ICA, or NNMF, or ICA.
% 
% If 'X' is a is a (time-by-stimuli) array of membrane potentials or 
%   spike-rates, a linear decomposition of 'DATA' is a pair of matrices 
%  'activations' (stimuli-by-nK) and 'response_waves' (time-by-nK) which
%   satisfies the equation: 
% 
%          X = (Y.response_waves * Y.activations') + Y.baseline
% 
% If 'X' is a data structure (such as the struct returned by utils.load),
%   the relevent data from 'X' is analysed (psth.wave or hekaData.PassData)
%   as analysis.linear is called from utils.load, this is designed to work
%   in conjuntion with that code (e.g. analyse spike data if -psth set). 
% 
% This analysis splits the input data into components which capture the
%   different features of the overall response. for PCA, these are the
%   directions of maximum variance, for ICA these are the directions of the
%   maximum departure from normality. 
% 
% -pca   : apply PCA decomposition to the responses. The response
%          components are the directions of maximum variance. 
% -nnmf  : apply NNMF decomposition to data. Good for spike-rates, 
%          see DOC NNMF. Response waves and activations are non-negative. 
% -ica   : apply ICA decomposition to data. 
% -iica  : apply ICA decomposition to transpose of data. 
% 
% Other Options: 
% -nK [6] : set number of components
% -rest   : enble/disable baseline computation. 
%           if X is an array, -rest [roi] enables the computation of the
%           baseline based on the median (or other quantile) of X in the
%           selected ROI. length(roi) = size(X,1). Example: 
%           > analysis.linear( spikerate, ..., '-rest', time < 0)
% 
% -quantile-rest [0.5]: quantile for resting calculation. median (0.5)
%                       unless analysing data.hekaData.PassData, in which
%                       case the default is 0.2 
% 
% Updated 28 August 2022 - Calvin Eiber <ceiber@ieee.org> 

% Parse Input arguments

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(v) varargin{find(named(v))+1};

if any(named('-nK')), 
  if any(named('-nK:')) % As string argument e.g. -nK:2
       nK = str2double(regexp(varargin{find(named('-nK'),1)},'\d+','match','once'));
  else nK = get_('-nK'); % e.g. -nK 2
  end
elseif nargin > 1 && isnumeric(varargin{1}), nK = varargin{1}; 
else   nK = 6; % default 6 components 
end

% Set Analysis mode
analysis_mode = 'PCA'; 
if any(named('nnmf')) || any(named('-nnmf')), analysis_mode = 'NNMF'; end
if any(named('ica')) || any(named('-ica')), analysis_mode = 'ICA'; end
if any(named('iica')) || any(named('-iica')), analysis_mode = 'IICA'; end

% For Baseline calculation:
q_rest = 0.5; % median

% Extract data to analyse
if isnumeric(data), X = data; % data was provided as input argument
  if any(named('-rest')), base_roi = get('-rest'); 
  else base_roi = false;
  end
elseif isstruct(data) % structure was provided to analyse 

  if any(named('-psth')), 
      X = data.psth.wave;
      base_roi = data.psth.time < 0; 
  elseif isfield(data,'hekaData')
      X = data.hekaData.PassData;     
      if isfield(data,'time'), base_roi = (data.time < 0); q_rest = 0.2; 
      else base_roi = false; 
      end
  else error('unclear how to interpret input structure')
  end
  if any(named('-rest')), base_roi = false; end % disable baseline
else error('unclear how to interpret input')
end

% Perform baseline computation if requested: 
if any(named('-quantile-rest')), q_rest = get_('-quantile-rest'); end
if any(base_roi), baseline = quantile(reshape(X(base_roi,:),[],1), q_rest);
else              baseline = 0; 
end

X = X - baseline; 

result = struct; 

if strcmp(analysis_mode,'NNMF')
    
    % if any(named('-symm'))
    % else 
    X_polarity = mode(reshape(sign(X(X~=0)),[],1));
   [score,coeff] = nnmf(max(X*X_polarity,0),nK);
    score = score * X_polarity; 
    coeff = coeff';
    
elseif strcmp(analysis_mode,'ICA')
    r = rica(X,nK,'Standardize',false); % max non-normality of SCORE
    coeff = r.TransformWeights; 
    score = X / coeff'; 
elseif strcmp(analysis_mode,'IICA')
    r = rica(X',nK,'Standardize',false); % max non-normality of COEFF
    score = r.TransformWeights; 
    coeff = (score \ X)';
else % default: PCA 
    [coeff,score] = pca(X,'Centered',false);
end

% Standard post-processing, rescale coeff to unit scale 

coeff = coeff(:,1:nK);  coe_scale = max(abs(coeff));
score = score(:,1:nK);

% flip the sign if both of the following are true: |min| > |max| and more
% negative than positive activations 
coe_sign = mode(sign(coeff)) .^ (abs(min(coeff)) > abs(max(coeff))); 
coe_sign(coe_sign == 0) = 1; 
coe_scale = coe_scale .* coe_sign; 

coeff = coeff * diag(1./coe_scale);
score = score * diag(coe_scale) ;



result.activations = coeff;
result.response_waves = score;
result.baseline = baseline;
result.response_scaleFactor = coe_scale;

return
