function [MI,MI_max] = mutualInformation(im1,im2,fast)
% Mutual Information calculation 
% (courtesy Brandon Munn, bmunn@sydney.edu.au)

persistent nD nI nP Y1 Y2

if isempty(nD) || nargin < 3 || ~fast

    nD = 10; 
    nI = 20; 
    nP = numel(im1); 
    [~,Y1] = hist(im1(:),nI);
    Y2 = quantile(im2(im2 > 0), (1:nD-1)/nD);
end

logsum = @(p) nansum(reshape(p,[],1));
im2 = arrayfun(@(y) sum(Y2 <= y), im2);

p_MI = zeros(nI,nD); 

for pp = 0:(nD-1), p_MI(:,pp+1) = hist(im1(im2 == pp),Y1); end

p_MI = p_MI / nP; 
p_M = sum(p_MI,1);
p_I = sum(p_MI,2);

MI = logsum( p_MI .* log2( p_MI ./ (p_I * p_M)) );
if nargout == 1, return, end

p_MM = diag(p_M); 
MI_max = logsum( p_MM .* log2( p_MM ./ (p_M'* p_M)) );
