% Beauchamp Cell 2020.m
%
% Beauchamp, M. S., Oswalt, D., Sun, P., Foster, B. L., Magnotti, J. F., Niketeghad, S., ... & Yoshor, D. (2020).
% Dynamic stimulation of visual cortex produces form vision in sighted and blind humans. Cell, 181(4), 774-783.

clear all;
c.efthr = 0.05; % what mag of electric field goes through the model, just a  speed thing
v.drawthr = 1;

PARALLEL = 0;
if PARALLEL
    numCores = feature('numcores');
    p = parpool(numCores);
end

%% define cortical and visual space
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-10, 80];
c.pixpermm = 8; c.e.radius = 0.25;
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-40, 40]; v.visfieldWidth= [-40,40]; v.pixperdeg = 12;
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);

%% figure 4
if 1% figure 4, YBN doing letters
    % folder CharacterDiscriminationYBN
    [all_v, xy] = Beauchamp_getDataFig4();
     letter(1).name = 'C';   letter(1).order = [5 1 3 8 12 20 18 ];
    letter(2).name = 'N';    letter(2).order = [3 8 12 20 11 5  9 17];
    letter(3).name = 'S';    letter(3).order =[1 3 8  12 11 10 9 13 17 18 19 20];
    letter(4).name = 'U';    letter(4).order = [8 12 20 19 14 9 5 ];
    n_electrodes = 20;
    trl.amp = 1.5;    trl.pw = 100*10^-6;    trl.dur =  50*10^-3;    trl.freq = 200;    trl.lag =  50*10^-3;
    lag = 100*10^-3; % time between each electrode
    trl.order = -1; 
else % 03-328, figure 6
    % folder CharacterDiscrimination03281
    trl.amp = 5;
    trl.pw = 100*10^-6;
    trl.dur =  100*10^-3;
    trl.freq = 120;
    lag=100*10^-3; % time between each electrode
end

% do the temporal pulse trains
for tt = 1:12
    trl.lag = trl.lag +lag;
    trl.simdur = 6;
    tp = p2p_c.define_temporalparameters();
    trl = p2p_c.define_trial(tp,trl);
    all_trl(tt) = p2p_c.convolve_model(tp, trl); % create response to the pulse
end

% locations of the phosphenes, taken from Fig. 2
figure(1); clf
for ee = 1:length(all_v.e)
    polarplot(all_v.e(ee).ang.*pi/180,  all_v.e(ee).ecc, 'ko'); hold on
    text(all_v.e(ee).ang.*pi/180, all_v.e(ee).ecc, num2str(ee));
    ax  = gca; ax.RLim = [0 13];
end

for ii = 1:n_electrodes
    disp([' generating maps for letter ', num2str(ii), 'out of ', num2str(n_electrodes)]);
    v.e = all_v.e(ii);
    c = p2p_c.define_electrodes(c, v);
    c = p2p_c.generate_ef(c);
    v = p2p_c.generate_corticalelectricalresponse(c, v); % create receptive field map for each electrode
    all_rfmaps(ii) = v; % generate the receptive field map for each electrode
end

for ex = 1; %:length(letter)
    figure(ex); clf
    p2p_c.plotretgrid((img./max(img(:)))*256, v, gray(256), ex, ['';]); drawnow; hold on
    % open the video
    filename = ['Beauchamp_Cell2020', letter(ex).name];  vid = VideoWriter([filename, '.avi']);
    vid.FrameRate = 30;    open(vid);
    for t = 0:(6*vid.FrameRate)-1  %for each video frame
        tmpt = 1+ceil((t/30)*tp.tSamp);
        img = zeros(size(squeeze(all_rfmaps(1).e.rfmap(:, :, 1))));
        for ii = 1:length(letter(ex).order) % add up the response of all of the electrodes at that moment in time
            trl = all_trl(ii);  % define the pulse train
            v = all_rfmaps(letter(ex).order(ii)); % define the electrode
            img = img + (mean(v.e.rfmap, 3)*trl.resp(tmpt));
        end
        p2p_c.plotretgrid(img*500, v, gray(256), ex); drawnow
        frame = getframe(gca);
        writeVideo(vid,frame);
    end
    close(vid);
end



% 
%         %% creating movies
% 
% 
% 
%         for t = 1:length(out(exp).val(1).resp_i)
%             % sim
%             img = zeros(size(squeeze(v.e(1).rfmap(:, :, 1))));
%             for ee=1:length(letter(exp).order)
%                 e = letter(exp).order(ee);
%                 img = img + max(v.e(e).rfmap,[],  3).*out(exp).val(end).resp_i(t); % just use the timing of an electrode simulated in the middle of the trial
%             end   p2p_c.plotretgrid(img*23, v, [], 1);
%          
%             t = text(-17, 15, 'simultaneous stimulation');
%             set(t, 'Color', [1 1 1])
%             frame = getframe(gca);
%             writeVideo(vid,frame);
%         end
%         for t = 1:length(out(exp).val(1).resp_i)
%             % seq
%             img = zeros(size(squeeze(v.e(1).rfmap(:, :, 1))));
%             for ee=1:length(letter(exp).order)
%                 e = letter(exp).order(ee);
%                 img = img + max(v.e(e).rfmap,[],  3).*out(exp).val(ee).resp_i(t); % just use the timing of an electrode simulated in the middle of the trial
%             end
%             p2p_c.plotretgrid(img*23, v, [], 1); hold on
%             t = text(-17, 15, 'sequential stimulation');
%             set(t, 'Color', [1 1 1])
%             frame = getframe(gca);
%             writeVideo(vid,frame);
%         end
%         close(vid);
%     end



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

    function [v,  xy] = Beauchamp_getDataFig4()
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

    function [v,  xy] = Beauchamp_getDataFig6()
    xy = [13.37	2.11; 11.44	2.14; 9.77	2.13; 8.02	2.19; 5.35	2.21];

    [theta,rho] = cart2pol(xy(:, 1),xy(:, 2));
    theta = theta * 180/pi;

    for ee = 1:size(xy, 1)
        v.e(ee).ang = theta(ee);
        v.e(ee).ecc = rho(ee);
    end
    end



