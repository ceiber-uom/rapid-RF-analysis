
function p = predictedResponse(data, varargin)
% prediction = predictedResponse( dat, [rdat], ... )
% 

named = @(s) strncmpi(s,varargin,numel(s));
get_ = @(v) varargin{find(named(v))+1};

% clear
% data = tools.load('20180703_Cell_1 #9','-dir','../HEKA Radon/MAT','-pca');

if nargin < 2
     f = figure; rdat = plots.plot_radon_IMG(data); close(f);
elseif any(named('-im')), rdat = get_('-im');
elseif isstruct(varargin{1}), rdat = varargin{1};
else f = figure; rdat = plots.plot_radon_IMG(data); close(f);
end


% if any(named('-im')), rdat = get_('-im');
% elseif any(named('-raw')), rdat = []; 
% else 
%     f = figure; rdat = plots.plot_radon_IMG(dat); delete(f);
% end

% compose visual stimuli in the rdat coordinate space

% stimulus string:
% {'S',[diam],[xy]}
% {'A',[d1 d2],[xy]}

stimuli = {'spot', [20:20:100 150:50:500], [0 0]};
stimuli(2,:) = {'annulus', (150:50:500), [0 0]};
stimuli{2,2}(2,:) = 100 * ones(size(stimuli{2,2}));

pattern_upsample = 4; 

upsample_range = linspace(rdat.range(1), rdat.range(2), ...
                          pattern_upsample*numel(rdat.range));
[gx,gy] = meshgrid(upsample_range); 

spot_ = @(c,r) double((gx-c(1)).^2 + (gy-c(2)).^2 <= (r).^2) ;

total_response = []; 

for s = 1:numel(stimuli)

    stim_type = stimuli{s,1};
    diameters = stimuli{s,2}';
    if size(stimuli,2) > 2 % explicit centerpoints
         center_xy = stimuli{s,3};
    else center_xy = [0 0];
    end

    nD = [size(diameters,1), size(center_xy,1)];

    for d = 1:max(nD)
        %% Make (antialiased) 

        diam = diameters( mod(d-1,nD(1))+1, :);
        c_xy = center_xy( mod(d-1,nD(2))+1, :);
        
        stim_pattern = spot_(c_xy, diam(1)/2);
        
        if upper(stim_type(1)) == 'A' % annulus 
            stim_pattern = stim_pattern - spot_(c_xy, diam(2)/2);
        end
        
        stim_pattern = imresize(stim_pattern,1/pattern_upsample); 
        % imagesc(rdat.range, rdat.range, stim_pattern)

        %% From stimulus image predict response wave
        
        % This is conceptual and will crash since I haven't got the
        % variable names correct

        nK = numel(rdat.images);
        c_img_kernal = rehsape(cat(3,rdat.images{:}), [], nK);
            
        % [1 x nPix] [nPix x nK] [nK x nT]
        c_score = reshape(stim_pattern,1,[]) * c_img_kernal;
        c_score = c_score + rdat.y_base; 

        total_response(:,end+1) = (c_score * rdat.profile)' ;
        
    end

end

total_response = total_response + dat.baseline; 

if nargout == 0 || any(named('-plot'))

    plot(dat.time, total_response)
    tidyPlotForIllustrator, xlim(dat.time([1 end]))
end

%%

error organise_output





return



function stimuli = parse_stim_pattern(get_,named)

