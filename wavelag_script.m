
clear

d = tools.load('20180703_Cell_1 #10','-dir','../HEKA Radon/MAT','-pca');
% 20180703 Mouse RF 1 1#3[on size ms]
% 20180703 Mouse RF 1#10[radon bars]

% r = plots.plot_radon_IMG(d); 

t = analysis.estimateWaveLag(d,'-plot'); 





