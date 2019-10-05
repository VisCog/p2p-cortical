% Fit_Winawer_data.m
%
% Finds parameters for v.drawthr and tp.scFac that fits phosphene area.



colorList =  [ 0.5 0 0 ;0.5 1 0.5; 1  0.8125 0;0 0.8750  1; 0  0 1.0000];

d = Winawer_getData(1:5);%% load the data

clear W
% pull out Winawer area and size, and set up 'site' structure
for ii=1:length(d)
    W{ii}.charge = d(ii).totalCharge.val;
    W{ii}.area = d(ii).poly_area.val;
    
    for tt = 1:length(d(ii).duration.val)
        site(ii).trl(tt).dur = d(ii).duration.val(tt); % duration in s
        site(ii).trl(tt).pw = d(ii).pulsewidth.val(tt) *10^-6; % pulse width in ms
        site(ii).trl(tt).freq = d(ii).frequency.val(tt); % -1 for a single pulse, NaN if not using a temporal model
        site(ii).trl(tt).amp = d(ii).current.val(tt)*1000; %should be in microAmps
        site(ii).trl(tt).draw.area = d(ii).area.val(tt);
        site(ii).trl(tt).draw.radius = d(ii).polycenterr.val(tt);
        site(ii).trl(tt).draw.brightness = d(ii).brightness.val(tt);
        site(ii).trl(tt).totalCharge = d(ii).totalCharge.val(tt);
    end
end



% Do the slow part and get v and c for each electrode to prepare for fitting
for ii=1:length(d)
    
    c.efthr = 0.01;
    v.drawthr = 0.05;
    % set up v and c for this electrode
    switch ii
        case 1
            v.e(1).ang = 19.8;      v.e(1).ecc = 26.6;  c.e(1).radius = 1.150 ; % in cm
            v.retinaSize = [30,30]; v.pixperdeg = 5;
            v.retinaCenter = [25,10];
        case 2
            v.e(2).ang = -166.4;    v.e(2).ecc = 9;     c.e(2).radius = 0.510 ;
            v.retinaSize = [20,20];v.pixperdeg = 8;
            v.retinaCenter = [-10,0];
        case 3
            v.e(3).ang = 142.2;     v.e(3).ecc = 5.12;  c.e(3).radius = 1.150 ;
            v.retinaSize = [8,8];v.pixperdeg = 8;
            v.retinaCenter = [-4,3];
        case 4 % central electrodes
            v.e(4).ang = 135;       v.e(4).ecc = 1.9;   c.e(4).radius = 1.150 ;
            
            v.retinaSize = [6,6];v.pixperdeg = 10;
            v.retinaCenter = [-1.25,1];
            
        case 5
            v.e(5).ang = 146.3;     v.e(5).ecc = 1;     c.e(5).radius = 1.150 ;
            v.retinaSize = [3.25,3.25];v.pixperdeg = 10;
            v.retinaCenter = [-.5,0];
            
    end
    c.cortexSize = [80,100];
    c.cortexCenter = [30,0];
    tp = p2p_c.define_temporalparameters();
    
    c = p2p_c.define_cortex(c);
    
    
    v = p2p_c.define_visualmap(v);
    [c, v] = p2p_c.generate_corticalmap(c, v);
    
    c = p2p_c.define_electrodes(c, v, ii);
    c = p2p_c.generate_ef(c, ii);
    % this is the slow part...
    v = p2p_c.generate_rfmap(c, v, ii);
    
    v_electrode(ii) = v;
    c_electrode(ii) = c;
end


% initial parameters
p.drawthr = .047;  %.05
p.scFac = 1;      % 1

getWinawerErr(p,v_electrode,c_electrode,tp,site,W);

bestP = fit('getWinawerErr',p,{'drawthr'},v_electrode,c_electrode,tp,site,W);





function d = Winawer_getData(sites)

pth = fullfile(ebsRootPath, 'data', 'ebs');
for ii = sites
    fname = sprintf('trial_data_site%d', ii);
    d(ii) = load(fullfile(pth, fname));
end
end