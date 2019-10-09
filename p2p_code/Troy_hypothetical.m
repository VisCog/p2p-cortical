% Troy_hypothetical.m
clear all; close all; clear mex

rng(1)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.efthr = 0.05; % what magnitude of electric field goes through the model, just a speed thing so one doesn't bother processing non-active regions of cortex
tp.scFac = 1;  % scaling of the strength of the pulse train
v.drawthr = .015; % what phosphene strength is visible
filename = 'Troyk_';
stim_method = 'sequential';
filename = [filename, stim_method];

%%
% set up the locations in terms of their positions on the cortical surface
i = 1;
for chip_x = 1
    for e_x = 1:10
        for e_y = 1:10
            c.e(i).radius = 1/1000; % scale is cm, if two small then the electrodes disappear due to undersampling of cortex
            c.e(i).hemi = 'lh';
            c.e(i).x = ((chip_x-1)*15) + ((e_x-1).*0.2);
            c.e(i).y = ((e_y-1).*0.2);
            i = i+1;
        end
    end
end

% define cortex & retina
c.cortexSize = [40,60]; % degrees top to bottom, degrees LR, divide by 2 to get the actual mm that are useful
c.pixpermm = 6; % default 6, for very small electrodes may need to be decreased
v.retinaSize = [40,40];v.pixperdeg = 5;  %visual field map size and samping
c = p2p_c.define_cortex(c); % build the cortical map
v = p2p_c.c2v_define_electrodes(c,v); % convert electrodes defined on the cortical surface to retinal co-ordinates
c = p2p_c.define_electrodes(c, v); % generates properties for each electrode
v = p2p_c.define_visualmap(v); % build the visual map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
c = p2p_c.generate_ef(c); % generate map of the electric field on cortical surface

tp = p2p_c.define_temporalparameters(); % define the temporal model


%% generate percepts for a pulse train
clist = parula(length(v.e));

for ii=1:length(v.e) % can work with parfor
    disp(sprintf('Electrode %d of %d',ii,length(v.e)));
    vtmp = p2p_c.generate_rfmap(c, v, ii); % generates the sum of weighted receptive fields activated by an electrod

    if strcmp(stim_method, 'sequential') trl.lag = (ii-1)*100*10^-3;% 50ms blank period between stimulation on each electrode
    else trl.lag = 100*10^-3; end % simultaneous
    trl.expname = 'Troyk_hypothetical'; 
    trl.e = ii;
    trl = p2p_c.define_trial(tp,trl); % define all the parameters of the trial, including the pulse train
    figure(1) ; 
    plot(trl.pt, 'Color', clist(ii,:)); hold on
    
    trl = p2p_c.generate_phosphene(vtmp, tp, trl); % multiply the rfmap by the response over time to create a percept over time
    % p2p_c.plotcortgrid(c.e(ii).ef * 256, c, gray(256),1,['subplot(1,2,1); title(''electric field'')']);
    % p2p_c.plotretgrid(trl.maxphos(:, :, 1)*100, v, gray(256), 1,['subplot(1,2,2); title(''phosphene'')']);
    savedtrial(ii) = trl;
end

save(filename, 'savedtrial'); % these are big files, so you may not want to save them, if you do, need to change
% your preferences/general/mat files to 7.3 or later

%% assuming simultaneous stimulation
% just adds together the max percept over time across all electrodes
figure(3); clf
img = savedtrial(1).maxphos(:, :, 1);
for ii=2:length(v.e)
    img = img+(savedtrial(ii).maxphos(:, :, 1)); hold on
    img(99:101, 99:101) = 255;
end
image(v.x, v.y, img); colormap(gray);


%% creating movie
close all
figure(1); clf
vid = VideoWriter([filename, '.avi']);
vid.FrameRate = 30;
open(vid);
itpl= round(linspace(1,length(savedtrial(end).pt), 30*savedtrial(end).trialdur));
for k = 1:length(itpl)
    img = zeros(size(squeeze(savedtrial(1).maxphos(:, :, 1))));
    for ii=1:length(savedtrial)
        if length(savedtrial(ii).resp)  >= itpl(k)
        img = img + squeeze(savedtrial(ii).maxphos(:, :, 1))./max(savedtrial(ii).resp).*savedtrial(ii).resp(itpl(k));
        end
    % this is a hack - maxphos = 1 * max response over time. So we divide
    % by that to get back to a phos with a max of 1, then multiply by the response over time.
    end
    img(99:101, 99:101) = 255;
    image(v.x, v.y, img*30); axis tight manual; axis off; colormap(gray);
    frame = getframe(gca);
    writeVideo(vid,frame);
end
close(vid)




