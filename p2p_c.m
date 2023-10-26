% p2p_c
%
%  holds all support functions for p2p cortex
% project.
%
% functions can be called from outside with 'p2p_c.<function name>'
% written IF and GMB
%
% 25/02/2023 moved into clean folder (IF)


classdef p2p_c
    methods(Static)
        %% definitions - temporal properties

        function trl = define_trial(tp, varargin)
            if nargin < 2;  trl = [];
            else; trl = varargin{1}; end
            if ~isfield(trl,'e'); trl.e = 1; end
            % durations are all in s
            if ~isfield(trl, 'dur');    trl.dur = 1000*10^-3;   end % duration in s of the electrical stimulation
            if ~isfield(trl, 'simdur');    trl.simdur = 3 ;   end
            % duration in s of the length of time being simulated, needs to be
            % longer to allow time for the neural response

            trl.t = 0:tp.dt:trl.dur-tp.dt;
            if ~isfield(trl, 'pw');     trl.pw = .1 * 10^-3;    end
            if ~isfield(trl, 'ip');     trl.ip = 0;             end % interphase delay
            if ~isfield(trl, 'lag');    trl.lag = 2*trl.pw;     end % delay before the pulse train begins
            if ~isfield(trl, 'order');  trl.order = 1;          end % 1 = cathodic first, -1  = anodic first
            if ~isfield(trl, 'freq');   trl.freq = 60;          end %NaN if not using a temporal model
            if ~isfield(trl, 'amp');    trl.amp = 100;          end % current amplitude in microAmps

            trl = p2p_c.generate_pt(trl, tp);
            trl.CperTrial = (trl.amp/1000) * trl.dur * trl.freq * trl.pw*10.^3; % charge per trial
            trl.CperPulse = trl.pw * trl.amp/1000; % charge per pulse
        end
        function trl = generate_pt(trl, tp)
            if isnan(trl.freq)
                % if all you are interested in is space then don't use a
                % temporal model at all, space and time are separable and
                % it's much faster
                trl.pt = 1;
            else
                if trl.ip ==0
                    on =  mod(trl.t,1/trl.freq) <=(trl.pw*2); % turn it on on
                    off = mod(trl.t-trl.pw,1/trl.freq) <=trl.pw & on;  % '& on' hack added by gmb;
                    tmp  = trl.amp.*(on-(2*off));
                else
                    on =  mod(trl.t,1/trl.freq) <trl.pw;
                    delay =  (trl.pw+trl.ip); % time difference between on and off
                    off = mod(trl.t-delay,1/trl.freq) < trl.pw;
                    tmp  = trl.amp.*(on-off);
                end

                lag = round(trl.lag/tp.dt); % delay before beginning the pulse train (frames)
                trl.pt= zeros(1, lag+length(tmp));
                trl.pt(lag+1:lag+length(tmp))=tmp;

                trl.t = 0:tp.dt:(trl.dur+trl.lag); % include the lag
                trl.t = trl.t(1:end-1);

            end
            if trl.dur<trl.simdur % usually we simulate a little longer than the trial, to allow for the response
                trl.pt((end+1):round((trl.simdur/tp.dt))) = 0;
                trl.t = 0:tp.dt:trl.simdur-tp.dt;
            end
        end

        function c = define_cortex(c)
            %% cortical magnification,
            % typical log z transformation parameters (based on early
            % Schwartz model
            if ~isfield(c, 'efthr'); c.efthr = 0.05; end % electric field values below this are assumed to be zero
            if ~isfield(c, 'animal') ;  c.animal = 'human'; end
            if strcmp(c.animal, 'human')
                c.k = 15; %imcale
                if ~isfield(c, 'a'); c.a = 0.5; end %fovea expansion for human, macaque is 0.3
                c.shift = c.k*log(c.a);
                if ~isfield(c, 'squish');   c.squish = 1;  end % some cortices are just a little rounder than others, no judgment
                if ~isfield(c, 'cortexHeight'); c.cortexHeight = [-40,40]; end %[height in mm of cortex, 0 is midline)
                if ~isfield(c, 'cortexLength'); c.cortexLength = [-5, 80]; end %[length in mm of cortex, 0 is fovea)
                if ~isfield(c, 'pixpermm');  c.pixpermm = 8; end    % choose the resolution to sample in mm.
            elseif strcmp(c.animal, 'macaque')
                c.k = 5; %scale
                c.a = 0.3; % values set by eyeballing Toottell data
                c.shift = c.k*log(c.a);
                if ~isfield(c, 'cortexHeight'); c.cortexHeight = [-20,20]; end %[height in mm of cortex, 0 is midline)
                if ~isfield(c, 'cortexLength'); c.cortexLength = [-5,30]; end %[length in mm of cortex, 0 is fovea)
                if ~isfield(c, 'pixpermm');  c.pixpermm = 8; end    % choose the resolution to sample in mm.
            elseif strcmp(c.animal, 'mouse')
                c.k = 1/40; % scale Garrett, 2014 FOR V1 how many mm of cortex represents 1 degree of visual field
            end
            %% receptive fields
            if ~isfield(c,'ar'); c.ar = 0.25; end % aspect ratio elongated rfs, Ringach 2002, J. Neurophysiology
            if ~isfield(c, 'rfmodel');   c.rfmodel = 'ringach';  end

            if strcmp(c.animal, 'human')
                if ~isfield(c, 'rfsizemodel')
                    c.rfsizemodel = 'keliris';  % Estimating average single-neuron visual receptive field sizes by fMRI. Keliris et al. PNAS 2019,
                end
                if strcmp(c.rfsizemodel, 'keliris') % using Keliris electrophysiology from supplementary table 1
                    c.slope = 0.08; % 0.05; % in terms of sigma of a Gaussian
                    c.intercept = 0.16; % 0.69;
                    c.min = 0;
                    c.delta = 2; % to do with on off receptive fields
                elseif strcmp(c.rfsizemodel, 'bosking') %  Saturation in Phosphene Size with Increasing Current Levels Delivered to Human Visual Cortex, Bosking et al. J Neurosci 2017
                    c.slope =   9
                    05; % Bosking data is in terms of diameter, so take the values from Figure 4 (slope = 0.2620 and intercept  = 0.1787) and divide by 2
                    c.intercept = 0.01;
                    c.min = 0;
                    c.delta = 2; % to do with on off receptive fields
                elseif strcmp(c.rfsizemodel, 'winawer')
                    c.slope = .1667;
                    c.min = 1.11;
                    c.intercept =  0.0721;
                    c.delta = 2; % to do with on off receptive fields
                end
            elseif strcmp(c.animal, 'macaque')
                c.slope =  0.06; % in terms of sigma
                c.intercept = 0.42;
                c.min = 0;
                c.delta = 2; % to do with on off receptive fields
                elseif strcmp(c.animal, 'macaque')
                c.slope =  0.06; % in terms of sigma
                c.intercept = 0.42;
                c.min = 0;
                c.delta = 2; % to do with on off receptive fields
            elseif strcmp(c.animal, 'mouse') %
                c.intercept = 20;  % Check this ezgi
            end

            %% ocular dominance columns
            if ~isfield(c, 'sig'); c.sig = .5;  end
            if ~isfield(c, 'onoff_ratio'); c.onoff_ratio  = 0.8; end % off cells contribute less to perception than on cells
            % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex. sig determines the distribution of OD values. Default is 5.
            % The larger sig, the more the distribution tends toward 0 and 1.

            if strcmp(c.animal, 'human')
                c.ODsize = 0.863; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 3; % 3mm creates the initial OD and orientation maps
            elseif strcmp(c.animal, 'macaque')
                c.ODsize = 0.531; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 1.85; % 3mm creates the initial OD and orientation maps
            elseif strcmp(c.animal, 'mouse')
                c.ODsize = NaN; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = NaN; % 3mm creates the initial OD and orientation maps
            end

            % define the size and resolution of cortical and visual space parameters
            c.gridColor = [1,1,0];
        end
        function [c] = generate_ef(c, varargin)
            % generates an electric field for each electrode, currently assumes they
            % are on the surface

            % max electric field is normalized to 1
            idx = 1:length(c.e);
            if ~isfield(c, 'emodel')
                c.emodel = 'Tehovnik';      end
            for ii=1:length(idx)
                R=sqrt((c.X-c.e(idx(ii)).x).^2+(c.Y-c.e(idx(ii)).y).^2);
                Rd = R-c.e(idx(ii)).radius; Rd(Rd<0) = 0;
                pt_ef = ones(size(c.X));
                if strcmp(c.emodel, 'Tehovnik')
                    if ~isfield(c, 'I_0');     c.I_0  = 1; end
                    if ~isfield(c, 'I_k');     c.I_k  = 6.75; end % controls rate of decay of current spread in cm
                    pt_ef=c.I_0./(1+c.I_k*Rd.^2);
                else
                    warning('electric field model not specified in code');
                    %    pt_ef(R>c.e(idx(ii)).radius)=2/pi*(asin(c.e(idx(ii)).radius./R(R>c.e(idx(ii)).radius)));
                end
                c.e(idx(ii)).ef = uint8(255.*pt_ef./max(pt_ef(:)));
            end
        end
        % generate functions
        function [c, v] = generate_corticalmap(c, v)
            % define cortex meshgrid
            c.x = linspace(min(c.cortexLength),max(c.cortexLength), (max(c.cortexLength)-min(c.cortexLength))*c.pixpermm);
            c.y = linspace(min(c.cortexHeight),max(c.cortexHeight), (max(c.cortexHeight)-min(c.cortexHeight))*c.pixpermm);
            [c.X,c.Y] = meshgrid(c.x,c.y);
            sz = size(c.X);

            %% Make the orientation and OD maps by bandpassing random noise

            % Rojer and Schwartz' method of bandpassing random noise:
            % Rojer, A.S. and E.L. Schwartz, Cat and monkey cortical columnar patterns
            %modeled by bandpass-filtered 2D white noise. Biol Cybern, 1990. 62(5): c. 381-91.
            %Make random noise: complex numbers where the angle is the orientation
            Z = exp(sqrt(-1)*rand(sz)*pi*2);

            % filter the noise to create initial columns
            freq = 1/c.ODsize; %cycles/mm (Try zero for big columns)
            filtPix = ceil(c.filtSz*c.pixpermm);
            [X,Y] = meshgrid(linspace(-c.filtSz/2,c.filtSz/2,filtPix),linspace(-c.filtSz/2,c.filtSz/2,filtPix));
            R = sqrt(X.^2+Y.^2);
            FILT = exp(-R.^2/c.sig.^2).*cos(2*pi*freq*R);  %Gabor

            %Convolve z with the filter
            W = conv2(Z,FILT,'same');
            c.ORmap = angle(W);
            WX = gradient(W);
            Gx = angle(WX);
            c.ODmap = normcdf(Gx*c.sig);


            %% Make the on and off  maps by bandpassing random noise
            % The idea that these vary smootly is based on
            % Najafian, S., Koch, E., Teh, K.L. et al. A theory of cortical map formation in the visual brain.
            % Nat Commun 13, 2303 (2022). https://doi.org/10.1038/s41467-022-29433-y
            % in this paper the maps are related to orientation and ocular
            % dominance, but because that doesn't matter for our
            % simulations we're just generating new maps

            % filter the noise to create initial columns
            freq = 1/c.ODsize *2; %cycles/mm, doubling the frequency
            filtPix = ceil(c.filtSz/2*c.pixpermm);
            [X,Y] = meshgrid(linspace(-c.filtSz/2,c.filtSz/2,filtPix),linspace(-c.filtSz/2,c.filtSz/2,filtPix));
            R = sqrt(X.^2+Y.^2);
            FILT = exp(-R.^2/c.sig.^2).*cos(2*pi*freq*R);  %Gabor

            %Convolve z with the filter
            W = conv2(Z,FILT,'same');
            u = (angle(W)/(pi)); % distance map
            tmp = zeros(size(u));
            tmp(u>0) = -log(u(u>0))/c.delta; % exponential fall of of d, as described by Mata & Ringach 2005
            tmp(u<0) = log(-u(u<0))/c.delta;
            c.DISTmap = tmp;
            WX = gradient(W);
            Gx = angle(WX);
            c.ONOFFmap = normcdf(Gx*c.sig);

            % create angle and eccentricity maps
            [c.v.X ,  c.v.Y] = p2p_c.c2v_real(c, c.X, c.Y);
            [c.v.ANG, c.v.ECC] = cart2pol(c.v.X,c.v.Y);
            c.v.ANG = c.v.ANG *180/pi;

            % create a mesh
            v.zAng = linspace(0,max(v.eccList),v.n)'*exp(sqrt(-1)*v.angList*pi/180);
            c.v.gridAng = p2p_c.v2c_cplx(c, v.zAng);

            v.zEcc = (v.eccList'*exp(sqrt(-1)*linspace(-90,90,v.n)*pi/180))';
            c.v.gridEcc = p2p_c.v2c_cplx(c, v.zEcc);
            c.RFsizemap = max(c.slope.* abs(c.v.ECC) + c.intercept, c.min); 

            if ~isfield(c, 'cropPix')
                c.cropPix  = c.v.ANG;
                c.cropPix(c.v.ECC>max([max(v.visfieldHeight) max(v.eccList)]))=NaN;
                if abs(min(c.cortexLength))<abs(max(c.cortexLength))
                    c.cropPix(c.X<0) = NaN;
                elseif abs(min(c.cortexLength))>abs(max(c.cortexLength))
                    c.cropPix(c.X>0) = NaN;
                end
            end
        end

        function v = define_visualmap(v)
            if ~isfield(v, 'visfieldHeight'); v.visfieldHeight = [-30 30]; end
            if ~isfield(v, 'visfieldWidth'); v.visfieldWidth = [-30 30]; end
            if ~isfield(v,'pixperdeg');     v.pixperdeg = 7;       end
            if ~isfield(v, 'drawthr');     v.drawthr = 1;    end
            v.x = linspace(v.visfieldWidth(1),v.visfieldWidth(2), (v.visfieldWidth(2)-v.visfieldWidth(1)).*v.pixperdeg);
            v.y = linspace(v.visfieldHeight(1),v.visfieldHeight(2), (v.visfieldHeight(2)-v.visfieldHeight(1)).*v.pixperdeg);
            [v.X,v.Y] = meshgrid(v.x, v.y);

            %Make the grid in retinal coordinates
            if ~isfield(v, 'angList');   v.angList = -90:45:90;    end
            if ~isfield(v, 'eccList');  v.eccList = [1 2 3 5 8 13 21 34];    end
            v.gridColor = [1 1 0];
            v.n = 201;
        end
        function [c, v] = define_electrodes(c, v)
            %  takes in the position of the electrode in visual
            %  co-ordinates and pop them onto the cortical surface
            idx = 1:length(v.e);
            if ~isfield(c, 'e') || ~isfield(c.e, 'radius')
                for ii=1:length(idx);        c.e(idx(ii)).radius = 500/1000;       end
            end
            if ~isfield(c.e, 'shape')
                for ii=1:length(idx);     c.e(idx(ii)).shape = 'round';    end
            end
            if ~isfield(v.e, 'ang')  % if putting in x, y co-ordinates rather than ang and ecc which is the default
                for ii = 1:length(idx)
                    [a, e]= cart2pol(v.e(idx(ii)).x, v.e(idx(ii)).y);
                    v.e(idx(ii)).ang = a*180/pi;
                    v.e(idx(ii)).ecc = e;
                end
            end
            for ii = 1:length(idx)
                [v.e(idx(ii)).x, v.e(idx(ii)).y] = pol2cart(v.e(idx(ii)).ang*pi/180, v.e(idx(ii)).ecc);
            end
            for ii = 1:length(idx)
                c.e(idx(ii)).area = pi*(c.e(idx(ii)).radius.^2);
                [c.e(idx(ii)).x, c.e(idx(ii)).y] = p2p_c.v2c_real(c,v.e(idx(ii)).x, v.e(idx(ii)).y);
            end
        end
        function v = c2v_define_electrodes(c,v)
            % If you've defined electrodes on cortex, this projects them
            % into visual space
            idx = 1:length(c.e);
            for ii = 1:length(idx)
                [v.e(idx(ii)).x, v.e(idx(ii)).y]= p2p_c.c2v_real(c,c.e(idx(ii)).x, c.e(idx(ii)).y);
                [v.e(idx(ii)).ang,v.e(idx(ii)).ecc] = cart2pol(v.e(idx(ii)).x,v.e(idx(ii)).y);
                v.e(idx(ii)).ang = v.e(idx(ii)).ang*180/pi;
            end
        end
        function [v, c] = generate_corticalelectricalresponse(c, v)
            if ~isfield(c, 'rfmodel')
                c.rftype = 'ringach';
            end
            % generates the sum of weighted receptive fields activated by an electrode
            % normalized so the max is 1
            idx = 1:length(c.e);
            for ii = 1:length(idx) % for each electrode
                if min(c.e(idx(ii)).x)-1<min(c.X(:)) || max(c.e(idx(ii)).x)+1>max(c.X(:)) ...
                        || min(c.e(idx(ii)).y)-1<min(c.Y(:))  ||  max(c.e(idx(ii)).y)+1>max(c.Y(:))
                    error('electrode is either outside or too close to the edge of the cortical sheet');
                end
                rfmap = zeros([size(v.X), 2]); % percept that includes a cortical model
                ct = 0;
                for pixNum = 1:length(c.X(:)) % for each cortical location
                    if (pixNum/400000 ==round(pixNum/400000))
                        per =round(100*pixNum/length(c.X(:)));
                        disp(['generating cortical electrical response ', num2str(per), '% complete']);
                    end
                    if ~isnan(c.cropPix(pixNum)) && abs(c.e(idx(ii)).ef(pixNum)) > c.efthr * 255
                        RF = p2p_c.generate_corticalcell(double(c.e(idx(ii)).ef(pixNum)), pixNum, c, v);
                        if ndims(RF)==4
                            RF = squeeze(RF(:, :, :, 1)); % don't worry about inhibition
                        end
                        if ~isempty(find(isnan(RF(:)), 1))
                            disp('wtf')
                        end
                        rfmap(:, :, 1)  =   rfmap(:, :, 1) + RF(:, :, 1);
                        rfmap(:, :, 2)  =   rfmap(:, :, 2) + RF(:, :, 2);
                        ct = ct+1;
                    end
                end
                if ct <c.e(ii).radius*10
                    if ct ==0;    disp('WARNING! No pixels passed ef threshold.');
                    else;    disp('WARNING! Very few pixels passed ef threshold.');
                    end
                    disp(' try checking the following:');
                    disp('lowering c.efthr or increase stimulation intensity');
                    disp('checking location of electrodes relative visual map');
                    disp('check the sampling resolution of cortex is not too low');
                end
                v.e(idx(ii)).rfmap = rfmap./max(abs(rfmap(:)));
            end
        end
        function [v, c] = generate_corticalvisualresponse(c, v)
            if ~isfield(c, 'rfmodel')
                c.rftype = 'ringach';
            end
            c.target.R = NaN(size(c.X));      c.target.R  =       c.target.R (:);
            % convolve the of  receptive fields with the visual stimulus
            for pixNum = 1:length(c.X(:)) % for each cortical location
                if (pixNum/100000 ==round(pixNum/100000))
                    per =round(100*pixNum/length(c.X(:)));
                    disp(['generating cortical visual response ', num2str(per), '% complete']);
                end
                if ~isnan(c.cropPix(pixNum))
                    RF = p2p_c.generate_corticalcell(1, pixNum, c, v);
                    if ~isempty(find(~isnan(RF(:)), 1))
                        on =  squeeze(RF(:, :, 1, 1)).*v.target.img; % on response to target
                        %off =  squeeze(RF(:, :, 1, 2)).*v.target.img; % suppressive response to target
                        c.target.R(pixNum) =sum(on(:));
                        error('this code no work yet')
                    end
                end
            end
            c.target.R = reshape(c.target.R, size(c.X));
        end

        function G = Gauss_2D(v,x0,y0,theta,sigma_x,sigma_y)
            % Generates oriented 2D Gaussian on meshgrid v.X,v.Y
            aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
            bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
            cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
            G = exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
        end

        function RF = generate_corticalcell(ef, pixNum, c, v)
            x0 = c.v.X(pixNum); % x center
            y0 = c.v.Y(pixNum); % y center
            od = c.ODmap(pixNum);
            theta = pi-c.ORmap(pixNum);  %orientation
            sigma_x = c.RFsizemap(pixNum)*c.ar; % minor axis sd
            sigma_y = c.RFsizemap(pixNum); % major axis sd

            % calculates a rf for a given location. Usually 3d (x, y and
            % eye) but for the ringach model it's 4d (x, y, eye,
            % on/off)
            if strcmp(c.rfmodel, 'scoreboard')
                % scoreboard version
                G = ef * exp(-( (v.X-x0).^2/(0.0001) + (v.Y-y0).^2/(.00001)));
                RF(:, :, 1) = G;
                RF(:, :, 2) = G;
            elseif strcmp(c.rfmodel, 'smirnakis')
                aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
                bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
                cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
                G = ef * exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
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
                A = sqrt((sum(tmp(:)>0.2)/v.pixperdeg.^2));
                d = c.DISTmap(pixNum) * A; % select a random d value

                % now create the real on and off fields, that are centered
                % on different locations in space and have variable
                % amplitudes
                x_off = x0 + (d/2).*cos(theta); y_off = y0 - (d/2).*sin(theta);
                x_on = x0 - (d/2).*cos(theta); y_on = y0 + (d/2).*sin(theta);

                wplus = c.ONOFFmap(pixNum); % scales the on subunit:  hplus_on and hminus_on
                wminus = 1-c.ONOFFmap(pixNum); % scales the off subunit: hplus_off and hminus_off

                % on subunit for bright dots
                hplus_on =  exp( - (aa*(v.X-x_on).^2 + 2*bb*(v.X-x_on).*(v.Y-y_on) + cc*(v.Y-y_on).^2));
                hplus_on = hplus_on./max(abs(hplus_on(:))); % response to bright dots, from the on subunit
                hplus_on = hplus_on * wplus;

                % off subunit, suppression with dark dots
                hplus_off = exp( - (aa*(v.X-x_off).^2 + 2*bb*(v.X-x_off).*(v.Y-y_off) + cc*(v.Y-y_off).^2));
                hplus_off = hplus_off./abs(max(hplus_off(:))); % response to dark dots, from the off subunit
                hplus_off = hplus_off *wminus * c.onoff_ratio;

                hminus_off = -0.4 * hplus_off; % inhibitory response to bright dots, from the off subunit
                hminus_on = -0.4 * hplus_on;% inhibitory response to dark  dots, from the on subunit

                % add the off component for bright and dark dots, and scale relative
                % amplitudes
                % bright = hplus_on + hminus_off; % response to bright dots
                % dark = hplus_off + hminus_on; % response to dark dots

                excitatory = hplus_on  - c.onoff_ratio*hplus_off; %
                inhibitory =hminus_on  - c.onoff_ratio*hminus_off ; % produces brightness where dark dots inhibiting

                RF(:, :, 1, 1)  =  od*ef*excitatory; % excitatory response
                RF(:, :, 2, 1)  =  (1-od)*ef*excitatory;

                RF(:, :, 1, 2)  =  od*ef*inhibitory;
                RF(:, :, 2, 2)  =  (1-od)*inhibitory;
            else
                error('c.rfmodel model not recognized')
            end
        end

        function [trl,v] = generate_phosphene(v, tp, trl)
            % finds the phosphene corresponding to the brightest moment in
            % time
            if isnan(trl.freq)
                trl.maxphos = v.e(trl.e).rfmap; trl.resp = 1;% the scaling due to current integration
            else
                % cunning hack to deal with the fact that the nonlinearity is not
                % spatiotemporally independent. We run the linear
                % spatiotemporally independent part of the model through the
                % convolve, then multiple it over space and pass that through
                % the nonlinearity.
                tmp = tp.model;
                tp.model = 'linear';
                trl = p2p_c.convolve_model(tp, trl);
                trl.maxphos = v.e(trl.e).rfmap.*max(trl.resp);
                tp.model = tmp; % restore the nonlinearity
                trl.maxphos = p2p_c.nonlinearity(tp, trl.maxphos);
                trl = p2p_c.convolve_model(tp, trl); % recalculate the response, using the nonlinearity;
            end

            % calculate the size of the image
            trl.sim_area = (1/v.pixperdeg.^2) * sum(trl.maxphos(:) > v.drawthr)/2; % calculated area of phosphene based on mean of left and right eyes
            if ~isempty(trl.maxphos)
                for i=1:2 % left and right eye
                    p = p2p_c.fit_ellipse_to_phosphene(trl.maxphos(:,:,i)>v.drawthr,v);
                    trl.ellipse(i).x = p.x0;                            trl.ellipse(i).y = p.y0;
                    trl.ellipse(i).sigma_x = p.sigma_x;     trl.ellipse(i).sigma_y = p.sigma_y;
                    trl.ellipse(i).theta = p.theta;
                end
                % what rule to use to translate phosphene image to brightness?
                beta = 6; % soft-max rule across pixels for both eyes
                trl.sim_brightness = ((1/v.pixperdeg.^2) * sum(trl.maxphos(:).^beta)^(1/beta));  % IF CHECK
            else
                trl.sim_brightness = [];
            end
        end
        function tp = define_temporalparameters(varargin)
            if nargin==0
                tp = [];
            else
                tp = varargin{1};
            end
            if ~isfield(tp, 'dt');       tp.dt = .001 * 10^-3; end % 001 time sampling in ms, should be no larger than 1/10 of tau1
            if ~isfield(tp, 'tau1');   tp.tau1 =0.0003; end % fixed based on Nowak and Bullier, 1998
            if ~isfield(tp, 'refrac');   tp.refrac = 100; end; %50 ; end % extent of attenuation for refractory period
            if ~isfield(tp, 'delta');   tp.delta = 0.001; end%0.001 ; end % decay parameter for refractory period

            if ~isfield(tp, 'tSamp');   tp.tSamp = 1000;end % subsampling to speed things up
            % fit based on Brindley, 1968, Tehovnk 2004 estimates 0.13-0.24 ms
            if ~isfield(tp, 'tau2');   tp.tau2 =  0.025; end%0.15;  end     % 24-33 from retina % IFCHANGE
            if ~isfield(tp, 'ncascades');  tp.ncascades = 3;   end% number of cascades in the slow filter of the temporal convolution
            if ~isfield(tp, 'gammaflag');   tp.gammaflag = 1;   end            %  include second stage game

            % leak out of charge accumulation
            %   tp.flag_cl=0; % 1 if you want to charge to leak back out of the system
            %   tp.tau2_cl = tp.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca

            % nonlinearity response parameters
             if ~isfield(tp, 'sc_in');   tp.sc_in =0.5663;   end    
            if ~isfield(tp, 'model');  tp.model = 'compression'; end
            if strcmp(tp.model, 'compression')
                if ~isfield(tp, 'power'); tp.power =  15.5901; end % chosen cos max brightness rating
                if ~isfield(tp, 'sc_out'); tp.sc_out = 10; end % fit using Winawer brightness data
            elseif strcmp(tp.model, 'sigmoid')
                disp('using sigmoid semisaturation constant')
                tp.asymptote = 2000;
                tp.e50 = 500; % electrical semisaturation constant
            elseif strcmp(tp.model, 'normcdf')
                disp('using normcdf semisaturation constant')
                tp.asymptote = 1500;
                tp.mean = 750;
                tp.sigma = 175;
            elseif strcmp(tp.model, 'weibull')
                disp('using weibull semisaturation constant')
                tp.asymptote = 1000;
                tp.thresh = 600;
                tp.beta = 3.5;
            end
        end

        %% psychophysics
        function v = generate_visualtarget(v)
            [x, y] = pol2cart(v.e.ang, v.e.ecc-v.target.offset);
            v.target.img = sqrt(((v.X-x).^2)+((v.Y-y).^2))<v.target.rad;
        end

        function [err, thresh] = loopall_find_threshold(tp,T)
            %
            % Runs the 'conv' model to get thresholds based on trials in the table 'T'.
            % returns the SSE and thresholds.  Table must contain fields holding the
            % following parameters for each trial:
            % pw      pulse width (sec)
            % dur     trial duration (sec)
            % freq    pulse frequency (Hz)
            % amp     amplitude at detection threshold
            if ~isfield(tp, 'nReps')
                tp.nReps = 12;
            end
            thresh = NaN(1,size(T,1));
            for i=1:size(T,1)
                % define trial parameters based on values in the table
                clear trl;  trl.pw = T.pw(i);   trl.amp = 1;    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 3; %sec
                trl = p2p_c.define_trial(tp,trl);

                tp.sc_in_ind = 1;     % separate 'scFac' for each experiment
                if isfield(tp,'experimentList')  % set tp_thresh accordingly for this trial
                    experimentNum = find(strcmp(T.experiment{i},tp.experimentList));
                    if ~isempty(experimentNum)  % set thresh_resp for this experiment.
                        tp.sc_in_ind = experimentNum;
                    end
                end
                thresh(i)= p2p_c.find_threshold(trl,tp);
            end

            err = nansum((thresh-T.amp').^2);
            %    disp(sprintf('tau1 = %g, tau2 = %g, power = %5.2f err= %5.4f  sc_in= %5.4f',  tp.tau1,tp.tau2,tp.power,err, tp.sc_in));
            disp(fprintf('mean = %g, sigma = %g,  err= %5.4f\n',  tp.mean,tp.sigma,err));
            if isfield(tp,'experimentList')
                for i = 1:length(tp.experimentList)
                    disp(fprintf('%10s: %g\n',tp.experimentList{i},tp.sc_in(i)));
                end
            end
        end
        function [err, thresh] = loop_find_threshold(tp,T)
            %
            % Runs the 'conv' model to get thresholds based on trials in the table 'T'.
            % returns the SSE and thresholds.  Table must contain fields holding the
            % following parameters for each trial:
            % pw      pulse width (sec)
            % dur     trial duration (sec)
            % freq    pulse frequency (Hz)
            % amp     amplitude at detection threshold
            if ~isfield(tp, 'nReps')
                tp.nReps = 12;
            end
            thresh = NaN(size(T,1), 1);
            for i=1:size(T,1)
                % define trial parameters based on values in the table
                clear trl;  trl.pw = T.pw(i);   trl.amp = 1;    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 3; %sec
                trl = p2p_c.define_trial(tp,trl);
                thresh(i)= p2p_c.find_threshold(trl,tp);
            end

            if strcmp('amp',T.Properties.VariableNames)
                err = nansum((thresh-T.amp').^2);
                disp(['err = ',  num2str(round(err, 3))]);
            else
                err = NaN;
            end

            if isfield(tp,'experimentList')
                for i = 1:length(tp.experimentList)
                    disp(sprintf('%10s: %g',tp.experimentList{i},tp.sc_in(i)));
                end
            end
        end
        function err = fit_brightness(tp, T)
            [loop_trl] = p2p_c.loop_convolve_model(tp,T);
            y_est = [loop_trl.maxresp]; y = [T.brightness];
            y_est = reshape(y_est, length(y), 1);
            y = reshape(y, length(y), 1);
            ind = ~isnan(y_est) & ~isnan(y);
            err = sum((y(ind)- y_est(ind)).^2);
            disp(sprintf('tau1 =%5.4f, tau2 =%5.4f, power =%5.4f, sc_in =%5.4f, sc_out =%5.4f, sse = %5.4f',  ...
                tp.tau1, tp.tau2, tp.power,tp.sc_in, tp.sc_out, err));
        end

        function [loop_trl] = loop_convolve_model(tp,T)
            %
            % Runs the 'conv' model to get thresholds based on trials in the table 'T'.
            % returns the SSE and thresholds.  Table must contain fields holding the
            % following parameters for each trial:
            % pw      pulse width (sec)
            % dur     trial duration (sec)
            % freq    pulse frequency (Hz)
            % amp     amplitude at detection threshold

            for i=1:size(T,1)
                % define trial parameters based on values in the table
                clear trl;  trl.pw = T.pw(i);   trl.amp = T.amp(i);    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 3; %sec
                trl = p2p_c.define_trial(tp,trl);

                % define impulse response
                if isfield(tp,'tSamp')
                    if tp.tSamp~=1% down-sample the time-vectors
                        t = trl.t(1:tp.tSamp:end);
                    end
                else
                    t = trl.t;
                end
                dt = t(2)-t(1);
                h = p2p_c.gamma(tp.ncascades,tp.tau2,t);            % Generate the n-cascade impulse response
                tid = find(cumsum(h)*dt>.999,1,'first'); % Shorten the filter if needed to speed up the code.
                if ~isempty(tid)
                    h = h(1:tid);
                else
                    sprintf('Warning: gamma hdr might not have a long enough time vector');
                end
                trl.imp_resp = h;  % close enough to use h
                loop_trl(i) = p2p_c.convolve_model(tp, trl);
            end
        end
        function amp = find_threshold(trl, tp)
            % Find amplitudes at threshold with the convolve model.
            % takes in trial, tp, and optional fitParams
            % finds and returns the trl.amp for which the max output of the
            % model for that trial, trial.resp, is equal to fitParams.thr

            if ~isfield(tp, 'nReps')
                tp.nReps = 12;
            end
            % first find the lowest 'hi' response
            hi = 1;
            resp = 0;
            while resp<tp.thresh_resp
                hi = hi*2;
                trl.amp = hi;
                trl = p2p_c.define_trial(tp,trl);
                trl = p2p_c.convolve_model(tp, trl);
                if tp.gammaflag
                    resp = max(trl.resp);
                elseif  tp.probsumflag
                    resp = trl.pd;
                end
            end
            lo = 0;
            % then do the binary search
            for i = 1:tp.nReps
                trl.amp  = (hi+lo)/2;
                trl = p2p_c.define_trial(tp,trl);
                trl = p2p_c.convolve_model(tp, trl);

                if max(trl.resp(:)) > tp.thresh_resp
                    hi = trl.amp;
                else
                    lo = trl.amp;
                end
            end
            amp = (hi+lo)/2;
        end
        function trl = convolve_model(tp, trl)
            % Implements 'finite_element' using the closed-form solution to
            % the respose to a pluse   Can be faster than 'finite_element'.
            % Assumes square pulse trains.
            %
            % tSamp is the temporal sub-sampling factor. Since tau2 is
            % relatively long, we can get away with a coarser temporal
            % sampling for the last convolution stage.  tSamp of 1000 works
            % well. Advise comparing to tSamp = 1 to check for innacuracy.
            % Also advise comparing 'convolve_model' to 'finite_element'
            % model which should be consiered the ground truth.
            %
            % Note: model only returns 'R3', 'spike' and 'resp' as output.
            % R1 and R2 (rectified R1) timecourses are not generated, so
            % the R2 of 'simpleleakyintegrator' has to obtained through the
            % 'finite_element' function.
            %
            % 'tt' is also returned, which is the temporally subsampled 't'
            % vector. Good for plotting 'spike' and 'resp'.

            % written GMB 6/17/2022

            % Assume the pulse train, pt, is a sequence of discrete jumps
            % in current. Find the 'events' where the pulse train, pt,
            % jumps up or down.
            %     Rconvtmp =  zeros(1,tp.ncascades+1); CHECK WITH GEOFF

            % Since spikes are sparse, manually convolve the 'spikes' with
            % the impulse response function at lowet temporal
            % resolution

            if isfield(tp,'tSamp')
                if tp.tSamp~=1% down-sample the time-vectors
                    t = trl.t(1:tp.tSamp:end);
                end
            else
                t = trl.t;
            end
            ptid = find(diff(trl.pt))+1;

            % R will hold the values of R1 at the event times.
            Rtmp = zeros(1,length(ptid));
            wasRising = 0;
            % Loop through the events, calculating R1 at the end of the event
            % and add impulse responses when R1 peaks and is after the refractory period.
            spikeId = logical(size(Rtmp));
            for i=1:(length(ptid)-1)
            %    tNow = trl.t(ptid(i+1));
                delta = trl.t(ptid(i+1))-trl.t(ptid(i));  % time since last 'event'
                % Closed form solution to leaky integrator that predicts
                % R(i+1) from R(i), delta and tau1:
                Rtmp(i+1) = trl.pt(ptid(i))*tp.tau1*(1-exp(-delta/tp.tau1)) + ...
                    Rtmp(i)*exp(-delta/tp.tau1);

                % Add a spike if:
                % (1) R1 is going down since last event
                % (2) R1 was going up before that, and
                % (3) we're past the refractory period since the last spike
                if Rtmp(i+1)<Rtmp(i) && wasRising
                    spikeId(i) = 1; % check spike id identical in both loops IF CHECK
                    wasRising = 0;  % no longer rising
                %    lastSpikeTime = trl.t(ptid(i+1));
                else
                    wasRising =1;
                end
            end

            R1= Rtmp*1000;
            R1(R1<0)=0;
            if isempty(R1) % no spikes at alll
                trl.resp = zeros(1, length(t));
                trl.resp_lin = trl.resp;
                trl.tt = t;  % for plotting
                trl.maxresp = 0;
                trl.spikeWhen= NaN;
                trl.imp_resp = NaN;
                trl.spikes = NaN;trl.spikes_norefrac = NaN;
            else
                spikeId (R1<=0) = 0; % sometimes non-spikes identified as spikes
                % pull out only spike events
                trl.spikes = R1(spikeId);

                ptid = ptid(spikeId);
                if tp.gammaflag
                    %   if ~isfield(trl, 'imp_resp')  % danger - if pre-computed,
                    %   it wont change if tau2 changes.
                    dt = t(2)-t(1);
                    h = p2p_c.gamma(tp.ncascades,tp.tau2,t);            % Generate the n-cascade impulse response
                    tid = find(cumsum(h)*dt>=.999,1,'first'); % Shorten the filter if needed to speed up the code. IF CHANGE
                    if ~isempty(tid)
                        h = h(1:tid);
                    else
                        disp('Warning: gamma hdr might not have a long enough time vector');
                    end
                    trl.imp_resp = h;  % close enough to use h
                    %    end
                    impFrames = 0:(length(trl.imp_resp)-1);
                    resp = zeros(1,length(t)+length(trl.imp_resp));        % zero stuff out

                    %reduction in spikes by inter-spike intervals:
                    interspike = [1,diff(trl.t(ptid))];
                    trl.spikes_norefrac = trl.spikes;
                    trl.spikes = trl.spikes.*(1-exp(-tp.refrac*(interspike+tp.delta)));
                    for i=1:length(trl.spikes)
                        id = find(t>trl.t(ptid(i)),1,'first');
                        resp(id+impFrames)  =   ...
                            resp(id+impFrames) + trl.spikes(i)*trl.imp_resp;
                    end
                    trl.resp_lin = resp;
                    resp = p2p_c.nonlinearity(tp, resp);
                    trl.maxresp = max(resp); % detection when maxresp goes above a threshold
                else
                    resp= trl.spikes;
                    trl.maxresp = max(resp);
                end
                % save the time-course of the response for output
                trl.resp = resp(1:length(t));
                trl.tt = t;  % for plotting
                trl.spikeWhen= ptid;
            end
        end

        %% utilities
        function out = chronaxie(p,pw)
            out = p.amp./(p.tau*(1-exp(-pw/p.tau)));
        end
        function y = gamma(n,k,t)
            %   y=gamma(n,k,t)
            %   returns a gamma function on vector t
            %   y=(t/k).^(n-1).*exp(-t/k)/(k*factorial(n-1));
            %   which is the result of an n stage leaky integrator.

            %   6/27/95 Written by G.M. Boynton at Stanford University
            %   4/19/09 Simplified it for Psychology 448/538 at U.W.
            %
            y = (t/k).^(n-1).*exp(-t/k)/(k*factorial(n-1));
            y(t<0) = 0;
        end
        function y = nonlinearity(tp,x)
            if ~isfield(tp, 'sc_in');    tp.sc_in = 1;
            else
                sc_in  =  tp.sc_in;
            end
            % some of our favorite static nonlinearities:
            switch tp.model
                case 'sigmoid'
                    y = sc_in .* x.^tp.power./(x.^tp.power + tp.sigma.^2);
                case 'normcdf'
                    y = normcdf(x, tp.mean, tp.sigma);
                    y(y<0) = 0;
                case 'weibull'
                    y = sc_in*p2p_c.weibull(tp,x);
                case 'power'
                    y = sc_in*x.^tp.power;
                case 'exp'
                    y = sc_in*x.^tp.k;
                case 'compression'
                    y = tp.sc_in.*(tp.power.*tanh((x*tp.sc_in)/tp.power));
                case 'linear'
                    y = tp.sc_in.*x;
            end
        end
        function [p] = weibull(params, x)
            % [p] = Weibull(params, x)
            %
            % The Weibull function based on this equation:
            %
            % k = (-log((1-e)/(1-g)))^(1/b)
            % f(x) = 1 - ((1-g) * exp(-(k*x/t).^b))
            %
            % Where g is performance expected at chance, e is performance level that
            % defines the threshold, b is the slope of the Weibull function, and t is
            % the threshold
            %
            % Inputs:
            %   params      A structure containing the parameters of the Weibull
            %               function:
            %       b       Slope
            %       t       Stimulus intensity threshold as defined by 'params.e'.
            %               When x = 'params.t', then y = 'params.e'
            %       g       Performance expected at chance, proportion
            %       e       Threshold performance, proportion
            %
            %   x           Intensity values of the stimuli
            %
            % Output:
            %   p           Output of the Weibull function as a function of the
            %               intensity values, x

            % Written by G.M. Boynton - 11/13/2007
            % Edited by Kelly Chang - February 13, 2017
            % Edited by Ione Fine - February 22, 2017

            if ~isfield(params, 'g')
                params.g = 0.5;
            end
            if ~isfield(params, 'e')
                params.e = (0.5)^(1/3);
            end
            k = (-log((1-params.e)/(1-params.g)))^(1/params.b);
            p = 1 - ((1-params.g) * exp(-(k*x/params.t).^params.b));
        end

        %% hypothetical array stuff
        function a = Array_Sim_Location(c, a)
            % simulates the left visual field (vx should be negative) and right cortex (cx should be positive)
            % takes in parameters
            % c - cortical surface
            % a - description of the array
            % p  - plotting parameters
            %
            % creates various arrays:
            % regular_cortex  - electrodes even on visual cortex
            % regular_visual field - evenly spaced phosphenes in visual apce
            % optimal - electrode spacing based on phosphene size

            if nargin <2
                a.arrayStyle = 'optimal';
            end
            if ~isfield(a, 'eccLim'); a.eccLim = [0 32]; end % range covered by the array
            if ~isfield(c, 'slope') ;  c.slope = .08; end % using Keliris electrophysiology from supplementary table 1
            if ~isfield(c,'intercept') ;  c.intercept = .16;    end
            if ~isfield(a, 'rf_sz'); a.rf_sz = (c.intercept+a.eccLim*c.slope); end %  rf sizes as a funtion of eccentricity

            if strcmp(a.arrayStyle, 'optimal')
                % Packs phosphenes in visual space that vary in size as a parametric  function of eccentricity.
                %  Phosphene centers are then projected into
                % Schwartz cortical space to show how spacing of electrodes are less
                % densely packed near the fovea.

                if ~isfield(a, 'spaceFac'); a.spaceFac = 1; end % separation of the electrodes, used for the optimal array

                % Pack phosphenes by generating optimal spacing along the horizontal meridian by adding up
                % phosphene sizes
                xi = 0;        si = 0;      i = 1;
                while(xi(end)+si(end)<a.eccLim(2))
                    si(i) = a.spaceFac*(xi(i)*c.slope+c.intercept);
                    xi(i+1) =(xi(i)+a.spaceFac*(si(i)+c.intercept)/2)/(1-c.slope/2);
                    i=i+1;
                end
                si(i) = a.spaceFac*(xi(i)*c.slope+c.intercept);

                % make concentric rings of phosphene at each spacing distance
                % a.vx a.vy are the positions of the phosphenes of the array in visual space
                ni = floor(2*pi*xi./si)';
                vx = 0;   vy = 0;
                for i=1:1:length(xi)
                    ang = linspace(0,2*pi,ni(i)+1)';
                    vx = [vx;xi(i).*cos(ang)];
                    vy = [vy;xi(i).*sin(ang)];
                end
                vx = vx-.0001;   % Hack to get the foveal phosphene inside the plotting range
                a.nelect  = length(vx); % calculate how many electrodes we get with this spacing

            elseif  strcmp(a.arrayStyle, 'regular_visualfield')
                if ~isfield(a, 'nelect'); a.nelect =1880; end % corresponds to a spacing of 1 with fov 32 using optimal array
                spacing = ceil(sqrt(a.nelect));
                if mod(spacing,2)==1
                    spacing = spacing +1;
                end
                tmp = linspace(-a.eccLim(2), a.eccLim(2), spacing);
                [vx, vy] = meshgrid(tmp, tmp); % match resolution of optimal array
                vx = vx(:); vy = vy(:);

            elseif  strcmp(a.arrayStyle, 'regular_cortex')
                if ~isfield(a, 'nelect'); a.nelect = 1880; end % corresponds to a spacing of 1 with fov 32 using optimal array

                [c_Length,~] = p2p_c.v2c_real(c, -a.eccLim(2),0);            % find axis limits
                [~,c_Height] = p2p_c.v2c_real(c,0,-a.eccLim(2));           % find axis limits
                c_Length = c_Length*1.2;     c_Height = c_Height*1.2;
                c_spacing = sqrt((c_Length*c_Height)/(a.nelect/2));

                c_Length = 0:c_spacing:c_Length;
                c_Height = 0:c_spacing:c_Height;
                c_Height = [-fliplr(c_Height(2:end)),c_Height];
                [cx, cy] = meshgrid(c_Length ,c_Height ); % create a grid on cortex

                cx = cx(:); cy = cy(:);
                ok = p2p_c.isValidCortex(c,cx,cy);
                [vx, vy] =  p2p_c.c2v_real(c,cx(ok), cy(ok));
                vx = cat(1, vx, -vx);vy = cat(1, vy, vy); % flip to the other side of the visual field
            else
                error('array.arrayStyle not recognized')
            end
            a.vx = vx(:); a.vy =vy(:);
            disp(a.arrayStyle)
            disp(['# electrodes = ',num2str(length(a.vx))]);
            a.actual_nelect = length(a.vx);
            a.res = c.pixpermm;
            a.Tfull = table(vx, vy);
            disp(['computed number of electrodes  = ' , num2str(length(vx))]);
        end
        function sz  = ecc2sz(a, ecc)
            % Support function for Array_Sim_Location
            % interpolate data to get phosphene size at locations x
            ecc = min(ecc, max(a.eccLim));
            sz = interp1(a.eccLim,a.rf_sz, ecc);
        end
        function Array_Sim_Plot(a, c, varargin)
            % Plots the array in visual space and cortical co-ordinates
            if nargin<3
                p.electrodeSize = 1; % p represents plotting parameters
            else
                p = varargin{1};
            end

            if ~isfield(p, 'electrodeSize');       p.electrodeSize =1; end
            if ~isfield(p, 'eccLim'); p.eccLim =[ 0, 32]; end %  [0,2.^ceil(log2(a.eccLim(2)))]; round to the nearest power of two
            if ~isfield(p, 'eccList'); p.eccList = [0,2.^(1:log2(p.eccLim(2)))]; end % eccentricity ring locations
            if ~isfield(p, 'angList'); p.angList = [90.1,135,180,225,270]*pi/180; end
            if ~isfield(p, 'ncol'); p.ncol = 256; end % number of colors in the colormap
            if ~isfield(p, 'cmap_scale'); p.cmap_scale = .25;end % how much of the colormap to use
            if ~isfield(p, 'symShape'); p.symShape = 'square'; end % shape of the symbols on the cortical map, 'square' or 'circle'
            if ~isfield(p, 'plotNums'); p.plotNums = 1:4; end

            % size vs eccentricity
            figure(p.plotNums(1));     set(gcf, 'Name' , "Phos Size vs. Ecc")
            ecc = linspace(p.eccLim(1),p.eccLim(2),p.ncol+1);
            sz = p2p_c.ecc2sz(a, ecc); % interpolate to find phosphene sizes at eccentricities defined by x
            idx = find(ecc<=a.eccLim(2));
            plot(ecc(idx), sz(idx),'-','Color',[.5,.5,.5],"LineWidth",2);
            set(gca,'YLim',[0,max(sz)*1.1]);    grid
            xlabel('Eccentricity (deg)');     ylabel('Phosphene size (deg)');

            % draw the phosphenes in visual coordinates, left visual field
            figure(p.plotNums(2)); clf; set(gcf, 'Name', 'left visual field electrodes'); hold on
            id = a.vx<=0;
            rad = sqrt(a.vx.^2+a.vy.^2);
            a.sz = c.slope*rad + c.intercept; % receptive field size of each electrode
            cmap = hsv(p.ncol);            % define colors for each electrode/phosphene
            col = 0.9*ones(length(id),3);  % gray
            scFac = p.cmap_scale*(p.ncol)/a.rf_sz(2); %scaling to get the right range in the colormap
            col(id,:) = cmap(min(ceil(a.sz(id)*scFac),p.ncol),:);

            ecc = repmat(p.eccList,p.ncol+1,1);
            ang = repmat(linspace(pi/2,3*pi/2,p.ncol+1)',1,length(p.eccList));
            gridx = ecc.*cos(ang);       gridy = ecc.*sin(ang);
            plot(gridx,gridy,'k-','Color',[.5,.5,.5]);
            ang = repmat(p.angList,p.ncol+1,1);
            ecc = repmat(exp(linspace(-20,log(max(p.eccList)),p.ncol+1))',1,length(p.angList));
            gridx = ecc.*cos(ang);    gridy = ecc.*sin(ang);
            cx = cos(linspace(-pi,pi,61));        cy = sin(linspace(-pi,pi,61));   % unit circle
            for i=1:length(a.vx)
                if sqrt(a.vx(i).^2+a.vy(i).^2) <a.eccLim(2)*1.4 
                patch(a.vx(i)+a.sz(i)/2*cx,a.vy(i)+a.sz(i)/2*cy,col(i,:));
                end
            end
            plot(gridx,gridy,'k-','Color',[.5,.5,.5]); axis equal;
            

            m=abs(p.eccLim(2))*1.1;  set(gca,'xLim',[-m,m]);   set(gca,'YLim',[-m,m]);   axis equal

            % Draw the corresponding 'electrodes' in cortical space, right
            % hemisphere
            figure(p.plotNums(3));  clf;   set (gcf, 'Name', "right hemi electrodes"); hold on
            ecc = repmat(p.eccList,p.ncol+1,1);
            ang = repmat(linspace(90.1,270,p.ncol+1)'*pi/180,1,length(p.eccList));
            x = ecc.*cos(ang);   y = ecc.*sin(ang);
            [gridx,gridy] = p2p_c.v2c_real(c,x,y);
            plot(gridx,gridy,'k-','Color',[.5,.5,.5]);
            ang = repmat(p.angList,p.ncol+1,1);
            ecc = repmat(exp(linspace(-20,log(max(p.eccList)),p.ncol+1))',1,length(p.angList));
            x = ecc.*cos(ang);  y = ecc.*sin(ang);
            [gridx,gridy] = p2p_c.v2c_real(c,x,y);

            id = a.vx<0;
            p.z = a.vx+sqrt(-1)*a.vy;
            [a.cx,a.cy] = p2p_c.v2c_real(c,a.vx,a.vy);
            if strcmp(p.symShape, 'square')
                shape_x = cos(-pi/4:pi/2:5*pi/4);   shape_y = sin(-pi/4:pi/2:5*pi/4);
            elseif strcmp(p.symShape, 'circle')
                shape_x = cos(linspace(-pi,pi,61));  shape_y = sin(linspace(-pi,pi,61));
            end
            for i=1:length(a.vx)
                if id(i);    patch(a.cx(i)+p.electrodeSize/2*shape_x, a.cy(i)+p.electrodeSize/2*shape_y,col(i,:));  end
            end
            plot(gridx,gridy,'k-','Color',[.5,.5,.5]);
            [limx,limy] = p2p_c.v2c_real(c,[-.01,-.01],max(p.eccList)*1.2*[1,-1]);            % find axis limits
            set(gca,'XLim',[-3,limx(1)]);       set(gca,'YLim',limy*1.2);            axis equal

            figure(p.plotNums(4));   set(gcf, 'Name', 'Colorbar');   clf
            maxSz = ceil(max(a.sz)); %deg
            img = repmat((1:(ceil(maxSz*scFac)))',1,2);
            image(1,linspace(0,maxSz,size(col,1)),img);
            colormap(cmap)
            set(gca,'YDir','normal');     set(gca,'XTick',[]);      ylabel('Phosphene size (deg)'); set(gca,'FontSize',18);
            axis equal;       axis tight
            set(gcf,'PaperPosition',[1,1,1,1]);
        end
        function rfmaps = Array_Sim_GenerateRFmaps(a, v, c, params, range)
            % Creates rfmaps for optimally spaced electrodes.
            % takes in:
            % a: (optional)  -  minimally has to define the arrayStyle and
            % the number of electrodes in total
            % v:  visual field parameters (optional)
            % c: cortical parameters (optional)
            % params, containing the folowing
            % range: (optional) -  which values in T will be
            % calculated, used for manual parallelization
            % flag_flip: (optional) - flip RFs to save time, can be  no
            % flip ([ 0 0 ], up down [1 0] (top to bottom), left right [0 1] (left to right), or up down
            % left right from upper left visual field quadrant [ 1 1]
            % overwrite = 0; % skip already existing file (0) or overwrite
            % them (1)
            %
            % Once you've run this you can use these rfs to make movies using
            % Array_Sim_Movies.m

            v = p2p_c.define_visualmap(v);
            c = p2p_c.define_cortex(c);
            [c, v] = p2p_c.generate_corticalmap(c, v);
            range = range(range<=height(a.T));

            rfmaps = NaN(size(v.X, 1), size(v.X, 2), length(range));
            for i = 1: length(range)
                e = range(i);
                disp(['electrode  ', num2str(e), ' within range ', num2str(min(range)) , '-', num2str(max(range))]);
                v.e.x = a.T.vx(e);
                v.e.y = a.T.vy(e);
                c = p2p_c.define_electrodes(c, v);
                c = p2p_c.generate_ef(c);
                v = p2p_c.generate_corticalelectricalresponse(c, v); % create receptive field map for each electrode
                rfmap = mean(v.e.rfmap, 3);
                if ~isempty(find(isnan(rfmap(:)), 1))
                    disp('NaNs,electrode likely outside the cortical sheet')
                else
                    if params.plot
                        p2p_c.plotretgrid(500*(rfmap+min(rfmap(:))), v, gray(256), 2);
                        figure(1); p2p_c.plotcortgrid(c.e.ef, c);
                        drawnow;
                        figure(2); p2p_c.plotretgrid(255*rfmap./max(rfmap(:)), v);
                        drawnow;
                    end

                    tmp = (rfmap+params.scale(1))*params.scale(2);
                    if max(tmp(:))>255
                        error('rfmap values too high for scaling to uint8 using current parameters');
                    elseif min(tmp(:))<0
                        error('rfmap values too low for scaling to uint8 using current parameters');
                    end
                    tmp = uint8(tmp);
                    rfmaps(:, :, i) = tmp;
                end
            end
        end
        function frame = Array_Sim_Movie(a, params, m)
            % Simulates what a movie would look like using a given array

            %% open input video file
            vid_in = VideoReader(m.filename_in);
            vid_dur  = vid_in.NumFrames;

            %% open output video file
            vid_out = VideoWriter(m.filename_out);
            vid_out.FrameRate = vid_in.FrameRate;
            open(vid_out);
            sz = size(a.rfmaps);

            for k = 1:vid_dur% for each frame
                frame = readFrame(vid_in);
                if k>=m.keepframes(1) && k<=m.keepframes(2)
                    disp(['simulating frame ', num2str(k), ' out of ', num2str(vid_dur)]);
                    if size(frame, 1)<size(frame, 2) % if the image isn't square
                        off = 1+ ( (size(frame, 2)-size(frame, 1))/2);
                        frame = frame(:, off-1:size(frame, 2)-off, :);
                    end
                    if ndims(frame)==3 % turn it grayscale
                        frame = mean(frame, 3);
                    end

                    in_movie=imresize(frame, [sz(1), sz(2)], "bilinear");
                    in_movie  = in_movie./255; % scale the movie between 0 -1

                    img = zeros(sz(1:2)); 

                    for e = 1:sz(3) % for each receptive field
                        orig_rfmap =  (double(a.rfmaps(:, :, e))./params.scale(2))  - params.scale(1); % goes between  -1 and 1
                        sz_rfmap = sqrt((sum(abs(orig_rfmap(:))))./(params.pixperdeg^2));
                        orig_rfmap = orig_rfmap/sz_rfmap;
                        if isempty(find(isnan(orig_rfmap(:)), 1))
                            rfmap = orig_rfmap;
                            amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
                            img = img + (rfmap.*amp);
                            if params.flag_flip(1)==1
                                rfmap = flipud(orig_rfmap);
                                amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
                                img = img + (rfmap.*amp);
                            end
                            if params.flag_flip(2)==1
                                rfmap = fliplr(orig_rfmap);
                                amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
                                img = img + (rfmap.*amp);
                            end
                            if sum(params.flag_flip)==2
                                rfmap = fliplr(flipud(orig_rfmap));
                                amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
                                img = img + (rfmap.*amp);
                            end
                        end
                    end
                    % weirdly, because we're not really simulating the temporal part of the
                    % model that relates electrical stimulation to brightness we don't include the nonlinearity
                    % we assume we want the brightness to be related to the image intensity
                    img =  img +  params.offset;
                    figure(2)
                    subplot(1,2,1)
                    imagesc(in_movie(m.crop:end-m.crop, m.crop:end-m.crop)); colormap(gray(256));
                    axis equal;    axis tight;    axis off
                    subplot(1,2,2)
                    image(img(m.crop:end-m.crop, m.crop:end-m.crop)); colormap(gray(256));
                    axis equal; axis off;  axis tight
                    drawnow
                    frame = getframe(gca);
                    writeVideo(vid_out,frame);
                end
            end
            close(vid_out);
        end


        %% transforms
        % transforms
        function c = define_cortical_size(c, v)
            vv = max(v.visfieldHeight);
            [ c.corticalHeight(2), ~] = p2p_c.v2c_real(c, 0, vv);
            vv = min(v.visfieldHeight);
            [c.corticalHeight(1), ~] = p2p_c.v2c_real(c, 0, vv);

            vv = max(v.visfieldWidth);
            [c.corticalWidth(2), ~] = p2p_c.v2c_real(c, vv, 0);
            vv =min(v.visfieldWidth);
            [c.corticalWidth(1), ~] = p2p_c.v2c_real(c, vv, 0);
        end
        function [vx, vy, ok] = c2v_real(c, cx, cy)
            % takes in real x, y numbers on the cortical grid and returns
            % real x, y numbers in visual space
            [z, ok] = p2p_c.c2v_cplx(c,cx + sqrt(-1)*cy);
            vx = real(z);
            vy = imag(z);
        end
        function w = v2c_cplx(c, z)
            % takes in imaginary numbers in visual space and returns imaginary
            % positions in the cortical grid

            lvf = real(z)<0; % find points in the left visual field
            z(lvf) = z(lvf)*exp(-sqrt(-1)*pi);  % rotate lvf points 180 degrees
            w = (c.k*log(z + c.a))-c.shift; % Use the Schwartz!
            w(~lvf) = w(~lvf)*exp(-sqrt(-1)*pi);  % rotate rvf (~lvf) points back 180 degrees

            % squish
            w = real(w)+c.squish*sqrt(-1)*imag(w);
        end
        function [cx, cy] = v2c_real(c, vx, vy)
            % takes in real x, y numbers in visual space and returns real
            % x, y positions on the cortical grid
            z = p2p_c.v2c_cplx(c,vx + sqrt(-1)*vy);
            cx = real(z);
            cy = imag(z);
        end
        function [w,ok] = c2v_cplx(c, z)
            % takes in imaginary numbers in cortical space and returns
            % imaginary positions in visual space. Not all points in
            % cortical space are valid, so invalid cortical points are
            % returned as NaNs in visual space.  'ok' is a logical where 0
            % is invalid and 1 is valid.

            % unsquish
            z = real(z)+ sqrt(-1)*imag(z)/c.squish;

            lvf = real(z)>0;  % find points in left visual field (right cortex)
            z(~lvf) = z(~lvf)*exp(-sqrt(-1)*pi); % rotate rvf points 180 degres
            w = exp((z+c.shift)/c.k)-c.a;  % Undo the Schwartz!
            w(lvf) = w(lvf)*exp(-sqrt(-1)*pi);  % rotate lvf points back 180 degrees

            % set the invalid points to NaNs
            ok = p2p_c.isValidCortex(c,real(z),imag(z));
            %ok = ones(size(w));  % gmb this needs to be removed
            w(~ok) = NaN;
        end
        function ok = isValidCortex(c,cx,cy)
            % Determines which points (cx,cy) are valid in cortical space
            % Points that fall outside the C-shaped map are invalid.

            % The boundary is determined by points (0,vyb) in visual
            % cortex.  For a given y value in cortical space, cy, the
            % corresponding vyb on the boundary is:
            vyb = c.a*tan(abs(cy)/c.k);

            % The corresponding x value on the boundary in cortex, cyb, is:
            cxb = c.k*log(sqrt(vyb.^2+c.a^2))-c.shift;

            % A point in cortical space is valid if cx is greater than the
            % corresponding x point on the boundary, which is when cx>cxb,
            % and cy falls within the asyptotic range of the map (+/-
            % k*pi/2)
            ok = abs(cx)>cxb & abs(cy)<c.k*pi/2;
        end
        %% plotting functions
        % plotting functions
        function logx2raw(base, precision)
            % logx2raw(base, precision)
            % Converts X-axis labels from log to raw values.
            % Inputs:
            %   base           Base of log transform (default: e)
            %   precision      Number of decimal places (default: 2)
            %
            % Example:
            % x = linspace(-3,0,11);
            % plot(log(x), log(x.^2));
            % logx2raw();
            % logy2raw(); % should be tolerant to multiple calls
            %
            % Note:
            % - See also: logy2raw.m

            % 11/17/96       gmb wrote it.
            % 6/6/96	     gmb added precision argument
            % 01/30/02       gmb updated it to use cell arrays, and to use original
            %                xtick values instead of converting labels. This way,
            %                multiple calls to this function doesn't keep converting
            %                the axis.
            % Edited by Kelly Chang - February 18, 2017

            %% Input Control

            if ~exist('base', 'var')
                base = exp(1);
            end

            if ~exist('precision', 'var')
                precision = 2;
            end

            %% Calculate Log x-axis Labels

            precision = sprintf('%%%2.1ff', precision*1.1);
            origXTick = get(gca, 'XTick'); % log x-axis labels (raw)
            newXTick = base.^(origXTick); % convert to raw
            newXLabel = arrayfun(@(x) sprintf(precision,x), newXTick, ...
                'UniformOutput', false); % write new x-axis labels
            set(gca, 'XTickLabel', newXLabel); % set x-axis labels of current graph
        end
        function logy2raw(base, precision)
            % logy2raw(base, precision)
            %
            % Converts Y-axis labels from log to raw values.
            %
            % Inputs:
            %   base           Base of log transform (default: e)
            %   precision      Number of decimal places (default: 2)
            %
            % Example:
            % x = linspace(-3,0,11);
            % plot(log(x), log(x.^2));
            % logx2raw();
            % logy2raw(); % should be tolerant to multiple calls
            %
            % Note:
            % - See also: logx2raw.m

            % 11/17/96       gmb wrote it.
            % 6/6/96	     gmb added precision argument
            % 01/30/02       gmb updated it to use cell arrays, and to use original
            %                xtick values instead of converting labels. This way,
            %                multiple calls to this function doesn't keep converting
            %                the axis.
            % Edited by Kelly Chang - February 18, 2017

            %% Input Control

            if ~exist('base', 'var')
                base = exp(1);
            end

            if ~exist('precision', 'var')
                precision = 2;
            end

            %% Calculate Log x-axis Labels

            precision = sprintf('%%%2.1ff', precision*1.1);
            origYTick = get(gca, 'YTick'); % log y-axis labels (raw)
            newYTick = base.^(origYTick); % convert to raw
            newYLabel = arrayfun(@(x) sprintf(precision,x), newYTick, ...
                'UniformOutput', false); % write new y-axis labels
            set(gca, 'YTickLabel', newYLabel); % set y-axis labels of current graph
        end
        function plotcortgrid(img, c,  varargin)
            % plotcortgrid(img, c)
            % plotcortgrid(img, c, cmap,figNum, evalstr)
            % takes as input:
            %   cortical image
            %   the structure c that defines the cortical surface
            % optional arguments:
            %   colormap, figure number and a string to evaluate
            %  (e.g. ''title('''corticalsurface''')' or 'subplot(1, 2,1)';


            if nargin<3 || isempty(varargin{1});  cmap = gray(256);   else; cmap = varargin{1}; end
            if nargin<4 || isempty(varargin{2});  figNum = 1;  else; figNum = varargin{2}; end
            if nargin<5 || isempty(varargin{3});  evalstr = '';      else;  evalstr = varargin{3}; end

            if isfield(c,'cropPix')
                img(isnan(c.cropPix)) = NaN;
                img= img+2;
                cmap = [0,0,0;cmap];
            end

            eval(evalstr); colormap(cmap);
            if ~isempty(img)
                image(c.x, c.y, img); hold on
            end
            xlabel('mm'); ylabel('mm')
            set(gca,'YDir','normal');
            if mean(c.cortexLength)<=0
                plot(c.v.gridAng, '-', 'Color', c.gridColor); hold on
                plot(c.v.gridEcc, '-', 'Color', c.gridColor);
            end
            if mean(c.cortexLength)>=0 %
                % plotting both hemispheres symmetrically
                plot(-c.v.gridAng, '-', 'Color', c.gridColor);
                plot(-c.v.gridEcc, '-', 'Color', c.gridColor);
            end

            axis equal;  axis tight
            set(gca,'XLim',[min(c.x(:)),max(c.x(:))]);
            set(gca,'YLim',[min(c.y(:)),max(c.y(:))]);
            drawnow;
        end
        function plotretgrid(img, v, varargin)
            % plotretgrid(img, c)
            % plotretgrid(img, c, cmap,figNum, evalstr)
            % takes as input:
            %   retinal image
            %   the structure v that defines the retinal surface
            % optional arguments:
            %   colormap, figure number and a string to evaluate
            %  (e.g. ''title('''corticalsurface''')' or 'subplot(1, 2,1)';

            if nargin<3 || isempty(varargin{1});  cmap = gray(256);   else; cmap = varargin{1}; end
            if nargin<4 || isempty(varargin{2});  figNum = 1;        else; figNum = varargin{2}; end
            if nargin<5 || isempty(varargin{3});  evalstr = '';      else; evalstr = varargin{3}; end

            figure(figNum); hold on
            eval(evalstr);
            image(v.x, v.y, img); hold on

            colormap(cmap); set(gca,'YDir','normal');

            plot(v.zAng,'-','Color', v.gridColor); plot(v.zEcc,'-','Color', v.gridColor);
            plot(-v.zAng,'-','Color', v.gridColor); plot(-v.zEcc,'-','Color', v.gridColor);

            axis equal;  axis tight
            xlabel('degrees'); ylabel('degrees')
            set(gca,'XLim',[min(v.x(:)),max(v.x(:))]);
            set(gca,'YLim',[min(v.y(:)),max(v.y(:))]);
            drawnow;
        end
        function p = fit_ellipse_to_phosphene(img,v)
            M00 = sum(sum(img));
            M10 = sum(sum(v.X.*img));           M01 = sum(sum(v.Y.*img));         M11 = sum(sum(v.X.*v.Y.*img));
            M20 = sum(sum(v.X.^2.*img));          M02 = sum(sum(v.Y.^2.*img));
            p.x0 = M10/M00;         p.y0 = M01/M00;
            mu20 = M20/M00 - p.x0^2;      mu02 = M02/M00 - p.y0^2;             mu11 = M11/M00 - p.x0*p.y0;
            a = (mu20+mu02)/2;         b = .5*sqrt(4*mu11^2+(mu20-mu02)^2);
            lambda_1 = a+b;      lambda_2 = a-b;
            p.theta = -.5*atan2(2*mu11,mu20-mu02);
            p.sigma_x = 2*sqrt(lambda_1);        p.sigma_y = 2*sqrt(lambda_2);
        end
        function fillSymbols(h,colList)
            if ~exist('h', 'var');     h = get(gca,'Children');    end
            for i=1:length(h)
                if ~exist('colList','var');       col = get(h(i),'Color');
                else
                    if iscell(colList); col = colList{i};
                    else;     col = colList(i,:);  end
                end
                set(h(i),'MarkerFaceColor',col);
            end
        end
        function draw_ellipse(trl, figNum, spstr, varargin)
            figure(figNum); hold on
            eval(spstr);
            if nargin >3
                eye = varargin{1};
            else
                eye = 1;
            end

            if nargin>4
                lineColor = varargin{2};
            else
                lineColor = 'g';
            end
            theta = linspace(-pi,pi,101);

            for e=1:length(eye)
                r = sqrt( (trl.ellipse(eye(e)).sigma_x*trl.ellipse(eye(e)).sigma_y)^2./ ...
                    (trl.ellipse(eye(e)).sigma_y^2*cos(theta).^2 + trl.ellipse(eye(e)).sigma_x^2*sin(theta).^2));
                x = trl.ellipse(eye(e)).x+r.*cos(theta-trl.ellipse(eye(e)).theta);
                y = trl.ellipse(eye(e)).y+r.*sin(theta-trl.ellipse(eye(e)).theta);
                plot(x,y,'-','LineWidth',1,'Color',lineColor);
            end
        end
        function  err = fit_phosphene(p,trl, v)
            % finds the standard deviation and mean of the best fitting
            % 2D Gaussian, assumes circularity.
            phos = mean(trl.maxphos, 3); phos = phos./max(phos(:));% average across the two eyes
            pred = p2p_c.Gauss(p, v);  pred = pred./max(pred(:));

            err = sum((phos(:)-pred(:)).^2);
            if p.sigma<0; err =err +  abs(p.sigma)*10.^6; end
        end
        function out = Gauss(p,v)
            out = exp(-((v.X-p.x).^2+(v.Y-p.y).^2)/(2*p.sigma^2));
        end
        function [err,predcx,predcy,cx,cy] = fitElectrodeGrid(p,vx,vy)
            % Gives a fit of a projected set of phosphenes onto a 6x4
            % electrode grid by first projecting the phosphenes using the
            % v2c mapping function and then moving and rotating the grid.

            p.shift = p.k*log(p.a);  % This really should be built in

            % Project phosphenes into cortex (including squishing)
            [cx,cy] = p2p_Beauchamp.v2c_real(p, vx, vy);

            % Creat a 6x4 array
            [x0,y0] = meshgrid(-1.5:1.5,-2.5:2.5);
            x0 =x0(:)';
            y0 =y0(:)';

            % Expand by dx
            x0 = x0*p.dx;
            y0 = y0*p.dx;

            % Rotate by ang and shift by (xc,yc)
            rot = [cos(p.ang) sin(p.ang);-sin(p.ang) cos(p.ang)];

            M = rot*[x0;y0] + repmat([p.xc;p.yc],1,24);
            predcx = M(1,:);
            predcy = M(2,:);

            % Compare predicted to projected
            err = sum((predcx-cx).^2 + (predcy-cy).^2);
        end
    end
end




