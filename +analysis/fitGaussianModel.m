

function gaussModel = fitGaussianModel(dat, varargin )
% gaussModel = fitGaussianModel(data, ... )
% 
% Fit a circular gaussian RF to the supplied flashing-bar data 
%  (expected in PCA format but should work for others). 
% 
% Each fitted gaussian has a center (x,y), radius, and an amplitude for
%   each PCA component, plus a baseline value (stimulus-independend).
% For example, if a DoG model (2 gaussians) is requested for data with 
%   4 PCA components, the resulting output has 
%   - 2 gaussian centers (4 parameters, x+y)
%   - 2 gaussian radii 
%   - 4 baseline (stimulus-independent) PCA component activations
%   - 8 amplitudes (one for each gaussian and each PCA component).
% 
% A correction for the finite bar width is applied: The bar width is
%   subtracted from the estimated radius. Bar width is assumed 1.5x delta-x
% 
% The output has the following fields: 
% .n_gaussians            : number of gaussians (nG)
% .fit_params   [nG x nP] : parameters of fit (unmodified)
% .center_xy    [nG x 2]  : gaussian center X and Y (typically in um)
% .gauss_radius [nG x 1]  : radius of each gaussian component
% .baseline     [nK x 2]  : baseline (stimulus-indepenent) activations of
%                           each PCA component 
% .amplitude    [nK x nG] : stimulus-dependent activations of each PCA
%                           component corresponding to each gaussian.
% .stats    : statistics structure with r2 and F-test results. The F-test
%             p-value indicates when the addition of an additional gaussian
%             component significantly reduces the residual variance, and is
%             only defined for nG > 1. 
% .resting  : waveform corresponding to the baseline (stimulus-independent)
%             PCA component activations.
% .kinetic  : waveform corresponding to the activation of each gaussian
%             component
% .predicted_RF : gaussian model receptive field activations for each
%                 stimulus (sinogram), should match the .y_all field output
%                 by utils.prepareRadon(). 
% 
% Options: 
% -nG [2]   : output a model with the selected number of components. 
%             if this option is not set, the output is an array of structs
%             corresponding to a 1-gaussian fit, a 2-gaussian fit, etc. 
% -nmax [2] : output models up to N fitted gaussians (default: 2 or DoG)
% -rmax [1] : allow gaussian radii up to X times the range of the input
%             data. Values of 0.0-2.0 are sensible. 
% -bw [1.5] : Apply bar width correction, where the bar width is X times
%             the difference in position for different bars. for no bar 
%             width correction, use ( ..., '-bw',0 )
% -no-plot  : Suppress generation of figure
% -no-wave  : Suppress computation of response kinetic profile for each
%              gaussian component
% -no-base  : Do not include resting level on kinetic profile plot
% -no-check : skip utils.prepareRadon on input data
% 
% v1.0 - 8 September 2022 - Calvin Eiber <c.eiber@ieee.org>

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(n) varargin{find(named(n))+1};

if any(named('-no-check')), d = dat; 
else d = utils.prepareRadon(dat, varargin{:}); 
end

%% Basics about stimulus

nK = size(d.y_all,2); % number of PCA components 

delta = d.x;
theta = deg2rad(d.ori);

nP = 2*nK+3; % number of parameters per gaussian 

% the arrangement of parameters within the gauss is as such:
% par(1:2) = gaussian XY position
% par(3)   = gaussian radius (circular)
% par(4:(nK+3)) = baseline (set to zero for higher components)
% par((nK+4):nP) = amplitude (for each PCA component)

%% Parse inputs

max_n_gaussians = 2;
if any(named('-nmax')), max_n_gaussians = get_('-nmax');
elseif any(named('-ng')), max_n_gaussians = get_('-ng'); 
end

r_max = 1.0; 
if any(named('-rmax')), r_max = get_('-rmax'); end

do_plot = ~any(named('-no-p')); 

bar_width = 1.5;
if any(named('-bw')), bar_width = get_('-bw'); end
bar_width = bar_width * mean(diff(unique(d.x)));

do_kinetic = isfield(d,'wave') && ~any(named('-no-w')); 

%% Set up functions for modeling response
c2x_ = @(xy) xy(1)*cos(theta) + xy(2)*sin(theta); % xy to delta given theta
gaussian_ = @(p,w) w(:,1) + w(:,2)*exp( -((delta-c2x_(p(1:2)))./(p(3))).^2 )';
% sinogram_ = @(y) reshape(y,[],nO);

LB = [delta(1)*[1 1] bar_width min(d.y_all) -range(d.y_all)];
UB = [delta(end)*[1 1 r_max] max(d.y_all) range(d.y_all)];
opts = optimoptions('fmincon','display','off');

%%

gaussModel = [];

y_meas = d.y_all;
y_model = zeros(size(y_meas));

for nG = 1:max_n_gaussians
    
    % how much unexplained variance per PCA component?
    residual = sum((y_meas-y_model).^2) ./ sum(y_meas.^2) ;
    initial_est = get_initial_guess(d, y_meas - y_model);

    % from initial guess, weight initial par by unexplained variance
    p0 = (residual) * initial_est / sum(residual);
    p0 = [p0(1:3) reshape(initial_est(:,4:5),1,[])];
    % initial guess of x,y,r and per-component weights

    c_ = @(p,n) p( (1:3) + (n-1)*nP );
    w_ = @(p,n) reshape( p( (4:nP)+(n-1)*nP ), [], 2); 

    parts_ = @(p) arrayfun(@(n) gaussian_(c_(p,n), w_(p,n)), 1:nG,'unif',0); 
    sumof_ = @(y) sum( cat(3,y{:}), 3)';
    v_ = @(x) reshape(x,[],1);

    gof = @(p) mean( (y_meas(:) - v_(sumof_(parts_(p)))).^2 );

    LB = ones(nG,1) * LB(1,:);
    UB = ones(nG,1) * UB(1,:);

    if ~isempty(gaussModel)

        p0 = [gaussModel(end).fit_params; p0]; %#ok<AGROW> 

        p0(2:end, 4:(nK+3)) = 0; % baseline only for first gaussian
        LB(2:end, 4:(nK+3)) = 0; 
        UB(2:end, 4:(nK+3)) = 0; 
    end

    p0 = v_(p0'); % convert to vector

    % find least-squares solution
    p1 = fmincon(gof, p0, [],[],[],[], v_(LB'), v_(UB'), [], opts);

    %% Assemble output

    % y_guess = sumof_(parts_(p0));
    y_model = sumof_(parts_(p1));

    weights = arrayfun(@(g) w_(p1,g), 1:nG,'unif',0);
    
    this = struct;
    this.n_gaussians = nG;
    this.fit_params = reshape(p1,[],nG)';
    [~,seq] = sort(this.fit_params(:,3),'ascend'); 

    this.center_xy = this.fit_params(seq,1:2);
    this.guass_radius = this.fit_params(seq,3) - bar_width;

    disp('TODO - here - confirm this bar width correction code is correct')
    disp('       look back at the code in the Radon paper masterfile')

    this.baseline = weights{1}(:,1);
    
    weights = cellfun(@(w) w(:,2), weights,'unif',0);
    this.amplitude = [weights{seq}]; 
    this.stats = []; 

    if do_kinetic
        this.resting = d.wave * this.baseline;
        this.kinetic = d.wave * this.amplitude;
        if isfield(dat,'resting_potential')
            this.resting = this.resting + dat.resting_potential;
        end
    end
    this.predicted_RF = y_model;

    if do_plot
        %% Make incremental fit plot
        
        clf
        subplot(2+do_kinetic,1,1), imagesc(y_meas');  title('data')
        % subplot(3+do_kinetic,1,2), imagesc(y_guess'); title('initial')
        subplot(2+do_kinetic,1,2), imagesc(y_model'); title('fitted model')
    
        if do_kinetic, subplot(3,1,3)
            plot(d.time, this.kinetic), hold on
            if ~any(named('-no-b'))
            plot(d.time, this.resting, 'Color', [.5 .5 .5]);
            end
            plot(d.time([1 end]),[0 0],'Color',[0 0 0 0.3]);
            tidyPlotForIllustrator, xlim(d.time([1 end]))

            for ss = 1:dat.nStimuli % add stim bars to the waves plot
              rectangle('Position',dat.stim_bar(ss,0.1), ... 
                        'FaceColor',[0 0 0 0.3], 'EdgeColor','none')
            end
        end

        h = get(gcf,'Children');
        set(h,'CLim',[-1 1]*max(abs([h.CLim]))) 
        suptitle(sprintf('%d gaussians', nG)); 
        pause(0.05)

    end

    if isempty(gaussModel), gaussModel = this;
    else gaussModel(end+1) = this; %#ok<AGROW> 
    end
end

%% Compute F-test stats for model

for gg = 1:nG

    stats = struct; 
    stats.r2 = corr( gaussModel(gg).predicted_RF(:), y_meas(:) ).^2; 

    if gg == 1
        stats.p = -1; 
        stats.f_test = [];
        gaussModel(gg).stats = stats; %#ok<AGROW> 
        continue, 
    end

    rm1 = y_meas - gaussModel(gg-1).predicted_RF; 
    rm2 = y_meas - gaussModel(gg).predicted_RF;

    [~,p,ci,s] = vartest2(rm1(:), rm2(:));

    stats.p = p;
    stats.f_test = s;
    stats.f_test.ci = ci;

    gaussModel(gg).stats = stats; %#ok<AGROW> 
end


if any(named('-ng')) % merge stats structures

    sngi = num2cell(1:nG); % stat # gaussians (for tracking)
    stats = [gaussModel.stats];
    [stats.n_gaussians] = deal(sngi{:});

    gaussModel = gaussModel(end);
    gaussModel.stats = stats;  
end

clear sngi stats p s ci nG this h weights p0 p1 rm1 rm2

%%
return



function guess = get_initial_guess(dat, y_all)

if nargin == 1, y_all = dat.y_all; end

nK = size(y_all,2);
orientations = unique(dat.ori); 
nO = numel(orientations);

guess = zeros(nK,5,nO); 

for kk = 1:nK
 for oid = 1:nO

  ori = (dat.ori == orientations(oid)); 
  
  x = dat.x(ori);
  y = y_all(ori,kk);

  yb = quantile(y,[.1 .95]);

  yi = cumsum(y-yb(1))/sum(y-yb(1));
  ix = interp1(yi,x,0.5);

  theta = deg2rad(orientations(oid));

  guess(kk,1:2,oid) = ix * [cos(theta) sin(theta)];
  guess(kk,3,oid) = sum(y > yb(2)/2) / 2.355 * mean(diff(x)); % approximately FWHM
  guess(kk,4:5,oid) = yb;

 end
end

guess = mean(guess,3);

