

% some of the best files:

% 20180703_Cell_1 #9
% 20180803_Cell_1 #1

clear
d = utils.load('20180703_Cell_1 #9','-dir','../HEKA Radon/MAT','-pca');

plots.standardFigure('Name','Standard PCA analysis'), clf
rdat = plots.plot_radon_IMG(d); 

% d = utils.prepareRadon(d, '-append'); 
% d = analysis.inverseRadon(d); 

plots.standardFigure('Name','Latency estimate'), clf
t = analysis.estimateWaveLag(d.response_waves(:,1), d.time, d.expoData,'-plot'); 

%% Select time points for display in "total RF" figure

nK = size(d.response_waves,2); 
% selection based on maximum wave amplitude (in complex sense) of each
% component
window = (d.time > -0.1 & d.time < 0.9*max(d.time))'; 
[~,tt_points] = arrayfun(@(k) max(abs(hilbert(d.response_waves(:,k))) ... 
                                      .* window), 1:nK ); 
tt_points = sort([t.zero_index t.index tt_points]); 

%% Generate plot of 'total RF' at each time-point

plots.standardFigure('Name','Total RF at timepoint'), clf

d.response_baseline = mean(d.response_waves(d.time<=0, :));
d.response_waves = d.response_waves - d.response_baseline;
subplot(3,1,1)
for k = 1:nK
    plot(d.time, 1e3*d.response_waves(:,k)), hold on   
end
for ss = 1:d.nStimuli,
  rectangle('Position',d.stim_bar(ss,0.1), ... 
            'FaceColor',[0 0 0 0.3], 'EdgeColor','none')
end

axis tight
try tidyPlotForIllustrator, end
set(gca,'XTick',unique(round(d.time(tt_points),2)),'XTickLabelRotation',-90)
set(gca,'Position',get(gca,'Position') + [0 2 0 -1]/50)


for tt = 1:numel(tt_points)

    total_rf_at_timepoint = 0 * rdat.image; 

    for k = 1:nK 
        total_rf_at_timepoint = total_rf_at_timepoint + rdat.images{k} * ...
                                     d.response_waves(tt_points(tt), k); 
    end

    subplot(3,4,tt + 4)
    imagesc(rdat.range,rdat.range, 1e3*total_rf_at_timepoint)
    axis image xy off
    title(sprintf('t = %0.2f', d.time(tt_points(tt))))
    set(gca,'Position',get(gca,'Position') + [-1 -1 2 2]/50)
end

h = flipud(get(gcf,'Children')); h(1) = []; 
set(h,'CLim',[-1 1] * max(abs([h.CLim])))

c = redbluecmap(11);
c = interp1((-5:5)', redbluecmap, linspace(-5,5,101)); 
colormap(c), 
ch = colorbar('southoutside'); 
ch.Position = [0.11 0.05 0.8 0.015];
xlabel(ch,'mV')

% In principle you could compute the radon transform of the measured data 
% at these time-points in particular to compare to this visualisation,
% which is based repeated rounds of baseline subtraction. I'm leaving this
% as an exercise to the next person to tackle this (which might very well
% be me, some time in the future)

