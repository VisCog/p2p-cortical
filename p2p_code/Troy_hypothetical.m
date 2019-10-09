% Troy_hypothetical.m
clear all; close all; clear mex
ampList = [1800]; % microAmps

rng(1)
c.efthr = 0.05; % what mag of electric field goes through the model, just a  speed thing
tp.scFac = 1;  % scaling of the strength of the pulse train
% v.drawthr = 0.15;
v.drawthr = .015;
filename = 'Troyk_';
stim_method = 'sequential';
filename = [filename, stim_method];
i = 1;
for chip_x = 1:4
    for e_x = 1:4
        for e_y = 1:4
            c.e(i).radius = 10/1000;
            c.e(i).hemi = 'lh';
            c.e(i).x = (chip_x*8) + ((e_x-1).*0.2);
            c.e(i).y = ((e_y-1).*0.2);
            i = i+1;
        end
    end
end


% define cortex & retina
c.cortexSize = [40,60]; c.pixpermm = 6; c.efthr = .05;
v.retinaSize = [40,40];v.pixperdeg = 5;  %10
c = p2p_c.define_cortex(c);
v = p2p_c.c2v_define_electrodes(c,v); % convert cortical electrodes to retina
c = p2p_c.define_electrodes(c, v);
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
c = p2p_c.generate_ef(c);

% create trials
trl.amp = 5000;
trl.freq = 50;
trl.pw = 2*10^(-4);
trl.dur = 1;
tp.model= 'normcdf';
tp = p2p_c.define_temporalparameters();

clist = parula(length(v.e));

%% loop
for ii=1:length(v.e)
    vtmp = p2p_c.generate_rfmap(c, v, ii);
    disp(sprintf('Electrode %d of %d',ii,length(v.e)));
    if strcmp(stim_method, 'sequential')
        trl.lag = (ii-1)*100*10^-3;% 50ms blank period between each stimulation
    else trl.lag = 100*10^-3; end
    trl.expname = 'Troyk_hypothetical';
    trl.e = ii;
    trl = p2p_c.define_trial(tp,trl);
    trl = p2p_c.generate_phosphene(vtmp, tp, trl);
    % p2p_c.plotcortgrid(c.e(ii).ef * 256, c, gray(256),1,['subplot(1,2,1); title(''electric field'')']);
    % p2p_c.plotretgrid(trl.maxphos(:, :, 1)*100, v, gray(256), 1,['subplot(1,2,2); title(''phosphene'')']);
    savedtrial(ii) = trl;
end

%% assuming simultaneous stimulation
figure(3); clf
img = savedtrial(1).maxphos(:, :, 1)*40;
for ii=2:length(v.e)
    img = img+(savedtrial(ii).maxphos(:, :, 1)*40); hold on
    img(99:101, 99:101) = 255;
end
image(v.x, v.y, img); colormap(gray);
%% creating movie
close all
figure(1); clf
vid = VideoWriter([filename, '.avi']);
vid.FrameRate = 30;
open(vid);
itpl= round(linspace(1,length(savedtrial(1).pt), 30*savedtrial(1).trialdur));
for k = 1:length(itpl)
    img = zeros(size(squeeze(savedtrial(ii).maxphos(:, :, 1))));
    for ii=1:length(v.e)
        img = img + squeeze(savedtrial(ii).maxphos(:, :, 1))./max(savedtrial(ii).resp).*savedtrial(ii).resp(itpl(k));
    end
    img(99:101, 99:101) = 255;
    image(v.x, v.y, img*30); axis tight manual; axis off; colormap(gray);
    frame = getframe(gca);
    writeVideo(vid,frame);
end
close(vid)
save(filename, 'savedtrial');



