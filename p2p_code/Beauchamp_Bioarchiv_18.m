% Beauchamp Bioarchiv 2018.m
%
% simulates Figure 2 from Beauchamp et al. BioRxiv 2018

clear all;
ampList = 800; % microAmps


PARALLEL = 1;

if PARALLEL
    numCores = feature('numcores');
    p = parpool(numCores);
end

[v, xy] = Beauchamp_getData();

c.efthr = 0.05; % what mag of electric field goes through the model, just a  speed thing
tp.scFac = 1;  % scaling of the strength of the pulse train
v.drawthr = .015;

c.cortexSize = [80,100]; c.pixpermm = 6;
v.retinaSize = [40,40]; v.pixperdeg = 5;


for ee = 1:length(v.e)
    c.e(ee).radius = 0.25;
end

% locations of the phosphenes, taken from Fig. 2
% figure(1); clf
% for ee = 1:length(v.e)
%     polarplot(v.e(ee).ang.*pi/180,  v.e(ee).ecc, 'ko'); hold on
%     text(v.e(ee).ang.*pi/180, v.e(ee).ecc, num2str(ee));     ax  = gca; ax.RLim = [0 13];
% end
%% build cortex and visual map
c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
c = p2p_c.define_electrodes(c, v);
tp = p2p_c.define_temporalparameters(tp);
c = p2p_c.generate_ef(c);
v = p2p_c.generate_rfmap(c, v);

%% order of stimulation based on figure 2
stim(1).order = [5 1 3 8 12 20 18 ];
stim(1).name = 'C';
stim(2).order = [3 8 12 20 11 5  9 17];
stim(2).name = 'N';
stim(3).order =[1 3 8  12 11 10 9 13 17 18 19 20];
stim(3).name = 'S';
stim(4).order = [8 12 20 19 14 9 5 ];
stim(4).name = 'U';

clist = parula(12);

parfor (exp = 1:length(stim), numCores * PARALLEL)
    electrode = generate_resp(stim(exp), tp, v)
    out(exp).val=electrode;
end

%% creating movies

for exp =  1:length(stim)
    figure(1); clf
    filename = ['beauchamp_', stim(exp).name];
    vid = VideoWriter([filename, '.avi']);
    vid.FrameRate = 30;
    open(vid);
    for t = 1:length(out(exp).val(1).resp_i)
        % sim
        img = zeros(size(squeeze(v.e(1).rfmap(:, :, 1))));
        for ee=1:length(stim(exp).order)
            e = stim(exp).order(ee);
            img = img + max(v.e(e).rfmap,[],  3).*out(exp).val(end).resp_i(t); % just use the timing of an electrode simulated in the middle of the trial
        end
        p2p_c.plotretgrid(img*23, v, [], 1);
        t = text(-17, 15, 'simultaneous stimulation');
        set(t, 'Color', [1 1 1])
        frame = getframe(gca);
        writeVideo(vid,frame);
    end
    for t = 1:length(out(exp).val(1).resp_i)
        % seq
        img = zeros(size(squeeze(v.e(1).rfmap(:, :, 1))));
        for ee=1:length(stim(exp).order)
             e = stim(exp).order(ee);
            img = img + max(v.e(e).rfmap,[],  3).*out(exp).val(ee).resp_i(t); % just use the timing of an electrode simulated in the middle of the trial
        end
        p2p_c.plotretgrid(img*23, v, [], 1); hold on
        t = text(-17, 15, 'sequential stimulation');
        set(t, 'Color', [1 1 1])
        frame = getframe(gca);
        writeVideo(vid,frame);
    end
    close(vid);
end



%% all the functions
function  electrode = generate_resp(stim, tp, v)

for ee=1:length(stim.order)
    trl = []; trl.expname = 'Beauchamp_BioRxiv';
    trl.e = stim.order(ee); trl.lag = (ee-1)*100*10^-3;% 50ms blank period between each stimulation
    trl = p2p_c.define_trial(tp,trl);
    trl = p2p_c.generate_phosphene(v, tp, trl);
    itpl= round(linspace(1,length(trl.resp), 30*trl.trialdur)); % interpolate to 30fps
    resp_i = trl.resp(itpl);
    electrode(ee).resp_i = resp_i;
end
end

function [v,  xy] = Beauchamp_getData()
xy = [5	8.678; 4.0805	7.842; 3.046	7.519; 2.3565	7.495; 6.4945	7.925; 5.2875	7.366; ...
    4.4255	7.049; 3.2185	6.6625; 7.931	6.653; 6.954	6.1595; 5.9195	5.951; 4.8275	5.741; ...
    8.9655	5.539; 7.701	5.3805; 6.8965	4.9505; 5.747	4.796; 10.23	3.916; 9.138	3.9935; ...
    7.701	3.7715; 6.5515	3.7315; 10.23	3.2265; 9.483	3.028; 8.3335	2.586; 6.954	2.309];

[theta,rho] = cart2pol(xy(:, 1),xy(:, 2));
theta = theta * 180/pi;

for ee = 1:size(xy, 1)
    v.e(ee).ang = theta(ee);
    v.e(ee).ecc = rho(ee);
end
end


