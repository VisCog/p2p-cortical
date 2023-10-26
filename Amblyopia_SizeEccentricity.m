
%% define cortical and visual space
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 80];
c.pixpermm = 12;
c.rfmodel= 'smirnakis';
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-10,10]; v.visfieldWidth= [0,60]; v.pixperdeg = 12;
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
% temporal parameters
tp = p2p_c.define_temporalparameters();

% set up plots
figure(1); hold on
p2p_c.plotretgrid(0, v, gray(64), 1, ['']);

trl.amp = [1000];
trl.dur =  200*10^-3;
trl.pw =  .1 * 10^-3;
trl.freq  = 200 * 10^-3;
trl = p2p_c.define_trial(tp, trl);
trl = p2p_c.convolve_model(tp, trl);
eccList = linspace(.1, 25, 5);
for ii=length(eccList)
    disp(sprintf('Electrode %d of %d',ii,length(eccList)));
    v.ecc = eccList(ii);  v.e.ang = 45; v.drawthr = 3;
    ind = find(c.v.ANG ==0 & c.v.ECC = v.ecc);
    ind = find(c.X>-1.5 & c.X<1.5 & c.v.>-1.5 & c.v.X<1.5
    RF= generate_corticalcells(c, v, 200);
    RF = RF./max(RF(:));
    p = p2p_c.fit_ellipse_to_phosphene(RF,v)
       p2p_c.plotretgrid(RF*256, v, gray(256), 2, ['';]);    drawnow; 
    %     trl.sim_radius= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
    %     trl.sim_diameter = 2 * trl.sim_radius; % we usually use radius/sigma, but Bosking uses diameter in this figure so we do too.
    %     trl.sim_brightness = max(trl.maxphos(:));
    %     sim_sizes(ii) =  trl.sim_diameter;
    %     sim_ecc(ii) = v.e.ecc;
end


    function sumRF = generate_corticalcells(c, v, nCells)
        sumRF = zeros(size(v.X));
        cInd = randperm(length(c.ODmap(:)));
       
        cInd  = cInd(1:nCells);       
        RFsize = max(c.slope .* v.ecc, c.min) + c.intercept;
        scatter =RFsize/2* randn(length(cInd), 2);
        for i = 1:nCells
            od = c.ODmap(cInd(i));

            theta = pi-c.ORmap(cInd(i));  %orientation
     
            sigma_x = RFsize* c.ar; % minor axis sd
            sigma_y = RFsize; % major axis sd
            x0 = v.ecc+scatter(i, 1); y0 = scatter(i, 2);

            % calculates a rf for a given location. Usually 3d (x, y and
            % eye) but for the ringach model it's 4d (x, y, eye,
            % on/off)
            if strcmp(c.rfmodel, 'scoreboard')
                % scoreboard version
                G = ef * exp(-( (v.X-x0).^2/(0.01) + (v.Y-y0).^2/(.01)));
                RF(:, :, 1) = G;
                RF(:, :, 2) = G;
            elseif strcmp(c.rfmodel, 'smirnakis')
                aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
                bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
                cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
                G =  exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
                G = G./sum(G(:));
                RF(:, :, 1)  = od*G;
                RF(:, :, 2)  = (1-od)*G;
            elseif strcmp(c.rfmodel, 'ringach')
                % creating more complex RFs based on Mata and Ringach, 2004 Neurophysiology paper
                % Spatial Overlap of ON and OFF Subregions and Its Relation to Response Modulation Ratio in Macaque Primary Visual Cortex
                aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
                bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
                cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);

                % calculate the area of of the on and off fields, which is
                % needed to normalize d. on and off fields use the same
                % oriented gaussian, only their central location and their
                % amplitudes differ
                tmp  = exp( -(aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
                A = sqrt(sum(tmp(:)>0.02)/v.pixperdeg.^2);
                d = c.DISTmap(cInd(i))/A; % select a random d value
                % now create the real on and off fields, that are centered
                % on different locations in space and have variable
                % amplitudes
                x_off = x0 + (d/2).*cos(theta); y_off = y0 - (d/2).*sin(theta);
                x_on = x0 - (d/2).*cos(theta); y_on = y0 + (d/2).*sin(theta);

                % on subunit for bright dots
                hplus_on =  exp( - (aa*(v.X-x_on).^2 + 2*bb*(v.X-x_on).*(v.Y-y_on) + cc*(v.Y-y_on).^2));
                hplus_on = hplus_on./sum(hplus_on(:)); % response to bright dots, from the on subunit
                hminus_on = 0.4 * hplus_on;% response to dark  dots, from the on subunit

                % off subunit with dark dots
                hplus_off = exp( - (aa*(v.X-x_off).^2 + 2*bb*(v.X-x_off).*(v.Y-y_off) + cc*(v.Y-y_off).^2));
                hplus_off = hplus_off./sum(hplus_off(:)); % response to dark dots, from the off subunit
                hminus_off = 0.4 * hplus_off; % response to bright dots, from the off subunit

                wplus = c.ONOFFmap(cInd(i)); % scales the on subunit:  hplus_on and hminus_on
                wminus = 1-c.ONOFFmap(cInd(i)); % scales the off subunit: hplus_off and hminus_off
                % add the off component for bright and dark dots, and scale relative
                % amplitudes
                bright = wplus*hplus_on + wminus*hminus_off; % brightness signal
                dark = wminus*hplus_off + wplus*hminus_on; % darkness signal

                RF(:, :, 1)  =   od*ef*(bright-c.onoff_ratio.*dark);
                RF(:, :, 2)  =  (1-od)*ef*(bright-c.onoff_ratio.*dark);
            else
                error('c.rfmodel model not recognized')
            end
            sumRF= sumRF+squeeze(sum(RF, 3));
        end
    end

