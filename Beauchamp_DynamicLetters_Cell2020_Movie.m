% Beauchamp_DynamicLetters_Cell 2020_Movie.m
%
% Replicates drawing task from Beauchamp dynamic stimulation paper
% Beauchamp, M. S., Oswalt, D., Sun, P., Foster, B. L., Magnotti, J. F., Niketeghad, S., ... & Yoshor, D. (2020).
% Dynamic stimulation of visual cortex produces form vision in sighted and blind humans. Cell, 181(4), 774-783.
% written IF & GMB
%
% loads in RF maps created by Beauchamp_DynamicLetters_Cell 2020_RF.m
% from a mat file Beauchamp_DynamicLetters_RF_Figure4.mat etc.
% 
% 25/02/2023 moved into clean folder (IF)
% 06/03/2023 Split RF generation and movie generation into separate files (IF) 

clear
debugflag = 0;
rng(11171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
fps = 30;
expList = { 'Figure 4-grid'}; %, 'Figure 6'}; % 'Figure 4'}; % grid is calculated location based on the array
Te = readtable("datasets\Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeLocations');
To = readtable("datasets\Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeOrder');

v.visfieldHeight = [-25, 25]; v.visfieldWidth= [0,25]; 
if debugflag
v.pixperdeg = 5;
c.pixpermm = 5; 

else
    v.pixperdeg = 12;
c.pixpermm = 12; 
end

v = p2p_c.define_visualmap(v);
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-60, 5];
c.e.radius = 0.25;
c = p2p_c.define_cortex(c);
[c, v] = p2p_c.generate_corticalmap(c, v);

flag = 'Simultaneous';
for ex = 1:length(expList)
    load(['datasets/Beauchamp_DynamicLetters_RF', expList{ex}]);
    % now for each letter, simulate phosphenes
    eid =strcmp(Te.experiment,expList{ex});
    Tloc = Te(eid, :);
    eid =strcmp(To.experiment,expList{ex});
    Torder = To(eid, :); % list of the letters/shapes simulated in this experiment

    for l =1:size(Torder, 1) % for each letter
        % open the video file for that letter
        filename = ['movies/Beauchamp_Cell2020_',expList{ex}, Torder.letter{l}, flag];
        vid = VideoWriter([filename, '.avi']);
        vid.FrameRate = fps;    open(vid);
        oList = str2num(Torder.order{l}); % order of stimulation for that letter
        img = zeros(size( saved(1).rfmap)); % zero out the frame
        trl.lag = .5;

        for o = 1:length(oList) % for each electrode,  as it is stimulated in order
            disp(['simulating electrode', num2str(o)]);
            disp(Tloc.electrode(oList(o)));
          trl.amp = Tloc.amp(o);    trl.pw = Tloc.pw(o);    trl.dur =  Tloc.dur(o);    trl.freq = Tloc.freq(o);
            if strcmp(flag, 'Sequential')
                trl.lag = 0.75 + ((o-1)*(Torder.pt_lag(l)+ Torder.pt_duration(l)));
            elseif strcmp(flag, 'Simultaneous')
                trl.lag = 0.75 ;
            end
            trl.simdur  = max(trl.lag) + .75;
            tp = p2p_c.define_temporalparameters();
            trl = p2p_c.define_trial(tp,trl);
            tp.model = 'linear'; % weird hack to deal with the spatiotemporal nonlinearity
            trl = p2p_c.convolve_model(tp, trl); % create response to the pulse
            t_ds = round(linspace(1, length(trl.resp), trl.simdur*fps)); % downsample response to video rate
            tp.model  = 'compression';
            for t = 1:length(t_ds)
                tmp =saved(oList(o)).rfmap*trl.resp(t_ds(t)); %  multiple space by time
                rfimg(o, t, :, :)   =single(p2p_c.nonlinearity(tp, tmp)); % pass it through the response nonlinearity
            end
        end
        for t = 1:length(t_ds)
            img = squeeze(sum(double(rfimg(:, t, :, :)), 1));
            p2p_c.plotretgrid(img*500, v, gray(256), ex); 
            drawnow; 
         set(gcf, 'Position',[1000         464        1247         874]);
            frame = getframe(gca);
            writeVideo(vid,frame);
        end
        close(vid);
    end % finished that letter
end % finished that experiment

