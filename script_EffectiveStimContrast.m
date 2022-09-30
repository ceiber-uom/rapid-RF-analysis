
clear
% n_reps = 32 random bootstrap centerings of RF 

% this script computes the effective contrast of a radon bar protocol and a
% checkerboard protocol

% this script also demonstrates the use of tools.artificalRadon 

X = linspace(-12,12,501); % image coordinates 
[gx,gy] = meshgrid(X);    % pixel grid

% equation of a 2D Circular Gaussian RF. 
% p[1] = center (x), p[2] = center (y), p[3] = radius
gauss_2d = @(x,y,p) exp( -( (x-p(1)).^2 + (y-p(2)).^2 ) ./ ...
                          (2*p(3)^2) ) / (2*pi*p(3)^2);

% DoG_RF - difference of gaussians 
DoG_RF = gauss_2d(gx,gy,[0 0 1]) - gauss_2d(gx,gy,[0 0 3]);

v_ = @(x) x(:); % column vector-ifier

dA = mean(diff(X)).^2; % size of a single pixel of DoG_RF

% Show RF map
figure(1), clf
imagesc(X,X,DoG_RF)
axis image xy

% title(sum(Z(:)) * dA)

if false 
    %% Look at response to a traditional (full-field) drifting grating
    disp('SpatFreq analysis') %#ok<UNRCH>     
    grating = @(x,f) cos(2*pi*f*x) + 1i*sin(2*pi*f*x); 
    
    sf = logspace(-2,0,101); % spatial frequency
    y_grat = arrayfun(@(f) v_(grating(gx(:),f))' * DoG_RF(:), sf) * (mean(diff(X)).^2);
    
    grating_max = max(abs(y_grat)); 
    absolute_max = sum(abs(DoG_RF(:))) * mean(diff(X)).^2;
    % semilogx(sf,abs(y_grat)), tidyPlotForIllustrator
end

w_value = logspace(-1,0.5,31);  % feature sizes to look at 
y_bars  = 0*w_value; 
y_check = 0*w_value;

clf, C = lines(7); W = @(i,v) (C(i,:)+v)/(1+v);

for ii = 1:length(w_value)

    rdat = tools.artificalRadon(DoG_RF,X,'width',w_value(ii),'no-figure');
    y_bars(ii) = max(abs(rdat.y)); % 2*.^2 / pi; 
    
    vx = linspace(rdat.x(1), rdat.x(end),32)'; % checkerboard X vector
    dx = mean(diff(vx))/2; % checkerboard pixel size
    
    [cx,cy] = meshgrid(vx); % coords of pixel centers for chekerboard
    checks = sign(randn([numel(vx).^2 numel(rdat.y)]));
    % do the same number of checkerboard patterns as radon patterns
    
    vc = zeros(size(rdat.y)); % RF drive per check pattern
    
    for uu = 1:numel(cx) 
        % For each pixel in cx, add the influence of that pixel across the
        % patterns
        im_x = (X >= cx(uu)-dx & X < cx(uu)+dx );
        im_y = (X >= cy(uu)-dx & X < cy(uu)+dx );        
        vc = vc + checks(uu,:)' * sum(v_(DoG_RF(im_x,im_y))) * dA; 
    end
    
    y_check(ii) = mean(abs(vc)); % average across patterns
    y_check_sd(ii) = std(abs(vc)); 
   
    %% Update the animation to show this bar and checkerboard
    
    clf
    imagesc(X,X,DoG_RF), axis image xy, hold on
    
    % XY coordinates of demo bar
    xy_bar = [-1 -1 1 1 -1; -1 1 1 -1 -1]' * [1.5 0; 0 21.5] * w_value(ii) * ...
           [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)] / 2 ;

    % XY coordinates of edges between pixels for checkerboard (lines)
    xy_check = [ v_([[1;1;nan]*vx' [1;-1;nan]*cy(1,:)]) ... 
                 v_([[1;-1;nan]*cy(1,:) [1;1;nan]*vx']) ];  

    % Show the demo bar against the simulated DoG RF in hazy blue
    plot(xy_bar(:,1),xy_bar(:,2),'-','LineWidth',1.1,'Color',W(1,0.5))
    
    % Show also the checkerboard in reddish 
    plot(xy_check(:,1),xy_check(:,2),'Color',W(2,0.5))
     
    axis tight, pause(0.05) % pause to show the figure
end

figure(2), clf
semilogx(w_value, y_bars), hold on
errorbar(w_value, y_check,y_check_sd,'CapSize',2)
try tidyPlotForIllustrator, end %#ok<TRYNC> 

xlabel('feature size')
ylabel('effective contrast')



