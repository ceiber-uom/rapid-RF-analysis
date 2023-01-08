

clear
% n_reps = 32 random bootstrap centerings of RF 

% this script computes the effective contrast of a radon bar protocol and a
% checkerboard protocol

% this script also demonstrates the use of tools.artificalRadon 

X = linspace(-120,120,501); % image coordinates 
[gx,gy] = meshgrid(X);    % pixel grid


% equation of a 2D Elliptical Gaussian RF. 
% p[1] = center (x), p[2] = center (y), p[3] = radius
gauss_2d = @(x,y,p) exp( -( ((x-p(1))/p(3)).^2 + ...
                            ((y-p(2))/p(4)).^2 + p(5) * ... 
                            ((x-p(1)).*(y-p(2))/(p(3)*p(4))) )); 

par = [15 -30 20 15 -0.3];
% ellipse_RF - elliptical gaussian
ellipse_RF = gauss_2d(gx,gy,par);

v_ = @(x) x(:); % column vector-ifier

dA = mean(diff(X)).^2; % size of a single pixel of DoG_RF

% Show RF map
figure(1), clf
imagesc(X,X,ellipse_RF)
axis image xy

hold on, plot(par(1),par(2),'r+')

% When the conic section is given in the general quadratic form
% A x2 + B xy + C y2 + Dx + Ey + F = 0

A = par(3)^-2;
B = par(5)./par(4)./par(3);
C = par(4)^-2;
F = -1;
D = 0;
E = 0;

% The following formula gives the eccentricity e if the conic section is
% not a parabola (which has eccentricity equal to 1), 
% not a degenerate hyperbola or degenerate ellipse, 
% and not an imaginary ellipse: 

eta = -sign(det( [A B/2 D/2 ; B/2 C E/2; D/2 E/2 F] )); 
ecc = sqrt( 2*sqrt((A-C)^2 + B^2) / ( eta*(A+C) + sqrt((A-C)^2+B^2) ));
disp(ecc)
title(sprintf('eccentricity = %0.3f', ecc))

%% generate synthetic radon data from elliptical RF

w_value = 8; % feature size to look at 
rdat = tools.artificalRadon(ellipse_RF,X,'width',w_value); % ,'no-figure');

%% Generate and view gauss fit to elliptical radon data

rdat.y_all = rdat.y - min(rdat.y) + 2.2;

figure(2), clf reset
gm = analysis.fitGaussianModel(rdat,'-ellipse','-nG',1);
disp(gm.gauss_eccentricity)

figure(3), plots.gaussianModel(rdat, gm)

% The fitted eccentricity doesn't look like it's exactly the same as the
% input eccentricity, not sure why this is but could have something to do
% with the bar width. 

%%

error TODO_generate_calibration_curve

% This is left as an exercise to the perfectly capable student
% 
% idea: 
% 
% in a loop, generate a bunch of different random RFs with eccentricity
% 0.01, 0.1, 0.2, ... 0.99. 
% make sure to get different radii and center positions of the RF
% 
% generate the simulated radon and fitted gauss model 
% 
% plot RF input eccentricity (determined from A, B, and C) and fitted
% eccentricity (gm.gauss_eccentricity) 

