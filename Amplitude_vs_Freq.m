% Amplitude vs. frequency

close all
clear all

debugflag = 0;

rng(3); %(1171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
% define pulse train
tp = p2p_c.define_temporalparameters(); % define the temporal model
trl.amp = 120; trl.freq = 200;
trl.pw = 2*10^(-4);   trl.dur= .8;

clear c
clear v
if debugflag
    v.pixperdeg = 5;  %visual field map size and samping
    c.pixpermm = 5; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
else
    v.pixperdeg = 14;  %visual field map size and samping
    c.pixpermm =14; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
end
scFac = 8;
c.cortexHeight = [-50, 50]; % degrees top to bottom, degrees LR,
c.cortexLength = [-40, 0];
v.visfieldHeight = [-2.5,2.5];
v.visfieldWidth= [0,5];
v.e.ecc  = 2.5;      v.e.ang = 0;
c = p2p_c.define_cortex(c); % define the properties of the cortical map
v = p2p_c.define_visualmap(v); % defines the visual map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
[c, v] = p2p_c.define_electrodes(c,v); % convert electrode locations from cortex to visual space
c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface
v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
figure(1); clf
freqlist = exp(linspace(log(1), log(2000), 11 ));
amplist = exp(linspace(log(1), log(400), 10));

ct = 1;
for f=1:length(freqlist)
    for a = 1:length(amplist)

        trl.amp = amplist(a);
        trl.freq = freqlist(f);
        trl = p2p_c.define_trial(tp,trl);
        % generate percepts for a pulse train
        tmp_trl = p2p_c.generate_phosphene(v, tp, trl);

        img = mean(tmp_trl.maxphos, 3);
        valbright(f,a) = tmp_trl.maxresp;
        valarea(f,a) =  tmp_trl.sim_area;
        img = (img*scFac) +125;
        ss = ['subplot(', num2str(length(freqlist)), ',', num2str(length(amplist)), ',', num2str(ct), ');'];

        p2p_c.plotretgrid(img,  v, gray(256), 1,ss);
        set(gca, 'FontSize', 6);
        titlestr = ['freq = ', num2str(freqlist(f)), ' amp = ', num2str(amplist(a))];
        title(titlestr)
        ct = ct + 1;
    end
end
figure(2); clf
surf(log(amplist), log(freqlist), valbright)
logx2raw;logy2raw;
xlabel('amp'); ylabel('freq');zlabel('brightness')
figure(3); clf
surf(log(amplist),log(freqlist),  valarea)

logx2raw;logy2raw;
xlabel('amp'); ylabel('freq');zlabel('area')


