function Array_Spacing_Movie()
%
% Replicates drawing task from Beauchamp dynamic stimulation paper
% Beauchamp, M. S., Oswalt, D., Sun, P., Foster, B. L., Magnotti, J. F., Niketeghad, S., ... & Yoshor, D. (2020).
% Dynamic stimulation of visual cortex produces form vision in sighted and blind humans. Cell, 181(4), 774-783.
% written IF & GMB
%
% This part
% 25/02/2023 moved into clean folder (IF)

clear
rng(11171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
fps = 30;


%% load image

in_filename = 'TeddySmiling.jpg';
in_img = imread(in_filename); in_img = mean(double(in_img), 3);
sz = min(size(in_img));
in_img = in_img(1:sz, 1:sz);
in_img_resize = imresize(in_img, [size(v.X, 1), size(v.X, 2)]);
out_filename = ['TeddySmiling', num2str(length(v.e))];


for o = 1:266% for each electrode,  as it is stimulated in order
    disp(['simulating electrode', num2str(o)]);
    load(['datasets/Optimal_Spacing/Optimal_Spacing_RFmaps_bosking_5e.mat');

    trl.simdur = 5;   trl.amp = Tloc.amp(o);    trl.pw = Tloc.pw(o);    trl.dur =  Tloc.dur(o);    trl.freq = Tloc.freq(o);
    if strcmp(flag, 'Sequential')
        trl.lag = 0.5 + ((o-1)*(Torder.pt_lag(l)+ Torder.pt_duration(l)));
    elseif strcmp(flag, 'Simultaneous')
        trl.lag = 0.5 ;
    end
    tp = p2p_c.define_temporalparameters();
    trl = p2p_c.define_trial(tp,trl);
    tp.model = 'linear'; % weird hack to deal with the spatiotemporal nonlinearity
    trl = p2p_c.convolve_model(tp, trl); % create response to the pulse
    t_ds = round(linspace(1, length(trl.resp), trl.simdur*fps)); % downsample response to video rate
    tp.model  = 'compression';
    for t = 1:length(t_ds)
        tmp =saved(oList(o)).rfmap*trl.resp(t_ds(t)); %  multiple space by time
        rfimg(o, t, :, :)   = p2p_c.nonlinearity(tp, tmp); % pass it through the response nonlinearity
    end
end
for t = 1:length(t_ds)
    img = squeeze(sum(rfimg(:, t, :, :), 1));
    p2p_c.plotretgrid(img*500, v, gray(256), ex); drawnow
    frame = getframe(gca);
    writeVideo(vid,frame);
end
close(vid);
end % finished that letter
end % finished that experiment
end
