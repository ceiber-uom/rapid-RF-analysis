
function hekaData = removeActionPotentials(hekaData) 

if ~isfield(hekaData,'rawPassData'), 
     hekaData.rawPassData = hekaData.PassData;
else hekaData.PassData = hekaData.rawPassData;
end

if isempty(hekaData.Spikes.spikeWaveforms)
    warning('No Spikes found!')
    return
end
    
spk_w = 32;

[coe,shape,lat] = pca(hekaData.Spikes.spikeWaveforms(1:spk_w,:),'Centered',false); 

coe = coe(:,1:3); 
shape = shape(:,1:3)';
lat = lat / sum(lat); %#ok<NASGU>

K = coe*shape;
K = K - K(:,1)*linspace(1,0,spk_w) - K(:,end)*linspace(0,1,spk_w);

n_spk = numel(hekaData.Spikes.timeIdx); 

for tt = 1:n_spk; 

    idx = hekaData.Spikes.timeIdx(tt) + (1:spk_w)-16; 
    pdx = hekaData.Spikes.passIdx(tt); 
    
    ok = (idx > 0 & idx < size(hekaData.PassData,1)); 
    idx = idx(ok); 
    
    hekaData.PassData(idx,pdx) = hekaData.PassData(idx,pdx) - K(tt,ok)'; 
%     
%     clf, set(gcf,'Color','w')
%     plot(time(1:spk_w), hekaData.rawPassData(idx,pdx)), hold on
%     plot(time(1:spk_w), hekaData.PassData(idx,pdx))
%     title(sprintf('%d/%d [p%d %0.1f s]', tt, n_spk, pdx, hekaData.Spikes.timeIdx(tt) / fs))
%     pause(0.1)
end

clear ok idx pdx tt n_spk spk_w K 