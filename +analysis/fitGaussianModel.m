

function estimates = fitGaussianModel(dat, varargin )

named = @(n) strncmpi(varargin,n,length(n));
get_ = @(n) varargin{find(named(n))+1};

if ~any(named('-no-check')), dat = utils.prepareRadon(dat, varargin{:}); end

%%

nK = size(dat.y_all,2);
orientations = unique(dat.ori); 
nO = numel(orientations);

nP = 2*nK+3; % number of parameters per gaussian 

%%

max_n_gaussians = 3;
if any(named('-nmax')), max_n_gaussians = get_('-nmax'); end

delta = dat.x;
theta = deg2rad(dat.ori);

c2x_ = @(xy) xy(1)*cos(theta) + xy(2)*sin(theta);
gaussian_ = @(p,w) w(:,1) + w(:,2)*exp( -((delta-c2x_(p(1:2)))./(p(3))).^2 )';
sinogram_ = @(y) reshape(y,[],nO);

r_max = 0.5; 
if any(named('-rmax')), r_max = get_('-rmax'); end

LB = [delta(1)*[1 1 0] min(dat.y_all) -range(dat.y_all)];
UB = [delta(end)*[1 1 r_max] max(dat.y_all) range(dat.y_all)];
opts = optimoptions('fmincon','display','off');

%%

iterative_fit = [];

y_meas = dat.y_all;
y_model = zeros(size(y_meas));

for nG = 1:max_n_gaussians
    
    residual = sum((y_meas-y_model).^2) ./ sum(y_meas.^2) ;
    initial_est = get_initial_guess(dat, y_meas - y_model);

    init = (residual) * initial_est / sum(residual);
    init = [init(1:3) reshape(initial_est(:,4:5),1,[])];


    c_ = @(p,n) p( (1:3) + (n-1)*nP );
    w_ = @(p,n) reshape( p( (4:nP)+(n-1)*nP ), [], 2); 

    parts_ = @(p) arrayfun(@(n) gaussian_(c_(p,n), w_(p,n)), 1:nG,'unif',0); 
    sumof_ = @(y) sum( cat(3,y{:}), 3)';
    v_ = @(x) reshape(x,[],1);

    gof = @(p) mean( (y_meas(:) - v_(sumof_(parts_(p)))).^2 );

    LB = ones(nG,1) * LB(1,:);
    UB = ones(nG,1) * UB(1,:);

    if ~isempty(iterative_fit)

        % init()

         error append_prev_fit
    end

    p_fit = fmincon(gof, init,  [],[],[],[], ...
                    v_(LB'), v_(UB'), [], opts);

    y_guess = sumof_(parts_(init));
    y_model = sumof_(parts_(p_fit));

    %%
    figure(1),clf
    subplot(3,1,1), imagesc(y_meas');
    subplot(3,1,2), imagesc(y_guess');
    subplot(3,1,3), imagesc(y_model');
    h = get(gcf,'Children');
    set(h,'CLim',[-1 1]*max(abs([h.CLim]))) 

    %%

    weights = arrayfun(@(g) w_(p_fit,g), 1:nG,'unif',0);

    this = struct;
    this.n_gaussians = nG;
    this.fit_params = v_(p_fit)';
    this.center_xy = [p_fit(1:nP:end); p_fit(2:nP:end)]';
    this.gaussian_radius = p_fit(2:nP:end);
    this.baseline = weights{1}(:,1);

    weights = cellfun(@(w) w(:,2), weights,'unif',0);
    this.amplitude = [weights{:}]; 

    if isempty(iterative_fit), iterative_fit = this;
    else iterative_fit(end+1) = this;
    end
    
end



















error todo



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

