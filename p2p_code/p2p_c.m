% p2p_c
%
%  holds all support functions for p2p cortex
% project.
%
% functions can be called from outside with 'p2p_c.<function name>'


classdef p2p_c
    methods(Static)
        % definitions
        function v = define_visualmap(v)
            
            if ~isfield(v, 'retinaSize')    v.retinaSize = [70,70]; end%  [height, width diameter in degrees]
            if ~isfield(v,'retinaCenter')    v.retinaCenter = [0,0];  end
            
            if ~isfield(v,'pixperdeg')
                v.pixperdeg = 7;
            end
            if ~isfield(v, 'drawthr')
                v.drawthr = 0.15;
            end
            v.x = linspace(.5/v.pixperdeg,v.retinaSize(2)-.5/v.pixperdeg,v.retinaSize(2)*v.pixperdeg)+v.retinaCenter(1) - v.retinaSize(2)/2;
            v.y = linspace(.5/v.pixperdeg,v.retinaSize(1)-.5/v.pixperdeg,v.retinaSize(1)*v.pixperdeg)+v.retinaCenter(2) - v.retinaSize(1)/2;
            [v.X,v.Y] = meshgrid(v.x, v.y);
            
            %Make the grid in retinal coordinates
            if ~isfield(v, 'angList')
                v.angList = -90:45:90;
            end
            if ~isfield(v, 'eccList')
                v.eccList = [1 2 3 5 8 13 21 34];
            end
            v.gridColor = [1 1 0];
            v.n = 201;
        end
        function c = define_cortex(c)
            %% cortical magnification,
            % typical log z transformation parameters (based on early
            % Schwartz model
            if ~isfield(c, 'animal')   c.animal = 'human'; end
            if strcmp(c.animal, 'human')
                c.k = 15; %scale
                c.a = 0.5; %fovea expansion for human, macaque is 0.3
            elseif strcmp(c.animal, 'macaque')
                c.k = 5; %scale
                c.a = 0.3; % values set by eyeballing Toottell data
            end
            %% receptive fields
            if ~isfield(c,'ar'); c.ar = 0.25; end
            if ~isfield(c, 'rfmodel');   c.rfmodel = 'smirnakis';  end
            if strcmp(c.rfmodel, 'smirnakis')
                if strcmp(c.animal, 'human')
                    c.slope =  0.05; % in terms of sigma
                    c.intercept = 0.69;
                    c.min = 0;
                elseif strcmp(c.animal, 'macaque')
                    c.slope =  0.06; % in terms of sigma
                    c.intercept = 0.42;
                    c.min = 0;
                end
            elseif strcmp(c.rfmodel, 'freeman')  % Freeman and Simoncelli, 2011
                c.slope = 0.1644; % receptive field diameter
                c.min = 1.21;
                c.intercept = 0;
            end
            %% ocular dominance columns
            if ~isfield(c, 'sig') c.sig = .5;  end
            % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex. sig determines the distribution of OD values. Default is 5.
            % The larger sig, the more the distribution tends toward 0 and 1.
            
            if strcmp(c.animal, 'human')
                c.ODsize = 0.863; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 3; % 3mm creates the initial OD and orientation maps
            else
                c.ODsize = 0.531; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 3; % 3mm creates the initial OD and orientation maps
            end
            
            % define the size and resolution of cortical and visual space parameters
            c.gridColor = [1,1,0];
            if ~isfield(c, 'cortexSize'); c.cortexSize = [80,100]; end %[height, width] Size of cortical maps (mm)
            if ~isfield(c, 'cortexCenter'); c.cortexCenter = [30,0]; end% center of electrode array (mm on cortex)
            if ~isfield(c, 'pixpermm');  c.pixpermm = 8; end    % choose the resolution to sample in mm.
            
        end
        function v = c2v_define_electrodes(c,v)
            idx = 1:length(c.e);
            for ii = 1:length(idx)
                v.e(idx(ii)).z = p2p_c.v2c(c,c.e(idx(ii)).x + sqrt(-1)*c.e(idx(ii)).y);
                v.e(idx(ii)).ang = angle(v.e(idx(ii)).z) * 180/pi;
                v.e(idx(ii)).ecc = abs(v.e(idx(ii)).z);
                v.e(idx(ii)).x = real(v.e(idx(ii)).z);
                v.e(idx(ii)).y = imag(v.e(idx(ii)).z);
            end
        end
        function c = define_electrodes(c, v)
            %  takes in the position of the electrode in visual co-ordinates
            idx = 1:length(v.e);
            
            if ~isfield(c.e, 'radius')
                for ii=1:length(idx)
                    c.e(idx(ii)).radius = 500/1000;
                end
            end
            if ~isfield(c.e, 'shape')
                for ii=1:length(idx)
                    c.e(idx(ii)).shape = 'round';
                end
            end
            
            for ii = 1:length(idx)
                c.e(idx(ii)).area = pi*(c.e(idx(ii)).radius.^2);
            end
            for ii=1:length(idx)
                if v.e(idx(ii)).ang>-90 && v.e(idx(ii)).ang <90
                    c.e(idx(ii)).hemi = 'lh';
                else
                    c.e(idx(ii)).hemi = 'rh';
                end
                % tranform electrodes
                
                for ii=1:length(idx)
                    if strcmp(c.e(idx(ii)).hemi, 'rh')
                        z = v.e(idx(ii)).ecc.*exp(sqrt(-1)*(v.e(idx(ii)).ang+180)*pi/180);
                    else
                        z = v.e(idx(ii)).ecc.*exp(sqrt(-1)*v.e(idx(ii)).ang*pi/180);
                    end
                    c.e(idx(ii)).z = p2p_c.c2v(c,z);
                    c.e(idx(ii)).x = real(c.e(idx(ii)).z);
                    c.e(idx(ii)).y = imag(c.e(idx(ii)).z); % turn into mm
                end
            end
        end
        
        function tp = define_temporalparameters(varargin)
            if nargin>0
                tp = varargin{1};
            end
            tp.dt = .001 * 10^-3; % time sampling in ms, should be no larger than 1/10 of tau1
            if ~isfield(tp, 'scFac'); tp.scFac = 1; end
            if ~isfield(tp, 'tau1')
                tp.tau1 =.012 * 10^-3;
            end% fit based on Brindley, 1968, Tehovnk 2004 estimates 0.13-0.24 ms
            tp.tau3 =  26.250* 10^-3; % 24-33 from retina
            
            % leak out of charge accumulation
            %   tp.flag_cl=0; % 1 if you want to charge to leak back out of the system
            %   tp.tau2_cl = tp.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca
            
            % nonlinearity parameters
            if ~isfield(tp, 'model');  tp.model = 'normcdf'; end
            if strcmp(tp.model, 'nanduri')
                disp('using nanduri semisaturation constant')
                tp.slope=3; %.3; % larger number = shallower slope
                tp.asymptote=14; %14;
                tp.shift=47; %47; % shifts curve along x-axis
            elseif strcmp(tp.model, 'sigmoid')
                disp('using sigmoid semisaturation constant')
                tp.asymptote = 2000;
                tp.e50 = 500; % electrical semisaturation constant
            elseif strcmp(tp.model, 'normcdf')
                disp('using normcdf semisaturation constant')
                tp.asymptote = 1500;
                tp.mean = 750;
                tp.sigma = 175;
            end
        end
        function trl = define_trial(tp,trl)
            %
            if ~isfield(trl,'expname'); trl.expname = 'generic'; end
            if ~isfield(trl,'e'); trl.e = 1; end
            
            if strcmp(trl.expname, 'Tehovnik')
                trl.dur = 100 * 10^-3; % duration in s
                trl.pw = .2 * 10^-3; % pulse width in s
                trl.order = -1; % 1 = cathodic first, -1  = anodic first
                trl.freq = 200; % -1 for a single pulse, NaN if not using a temporal model
                trl.amp = 50; % current amplitude in microAmps
            elseif strcmp(trl.expname, 'Bosking')
                trl.dur = 200*10^-3; % duration in s
                trl.t = 0:tp.dt:trl.dur-tp.dt;
                trl.pw = .1 * 10^-3; % pulse width in s
                trl.order = 1; % 1 = cathodic first, -1  = anodic first
                trl.freq = 200; % -1 for a single pulse, NaN if not using a temporal model
            elseif strcmp(trl.expname,'Evans')
                trl.order = -1;
                if ~isfield(trl, 'pw');     trl.pw = .25/1000;     end
                if ~isfield(trl, 'dur');    trl.dur = .5;   end% duration in ms
                if ~isfield(trl, 'freq');   trl.freq = 50;          end %NaN if not using a temporal model
                trl.t = 0:tp.dt:trl.dur-tp.dt;
            elseif strcmp(trl.expname, 'Beauchamp_BioRxiv')
                trl.order = -1;
                if ~isfield(trl, 'amp'); trl.amp = 1000; end
                if ~isfield(trl, 'pw');     trl.pw = .1*10^-3;     end
                if ~isfield(trl, 'dur');    trl.dur =  50*10^-3;   end% duration in ms
                if ~isfield(trl, 'freq');   trl.freq = 200;        end %NaN if not using a temporal model
                if ~isfield(trl, 'trialdur'); trl.trialdur = 1.7;  end % allows filler at end of trial
                trl.t = 0:tp.dt:trl.dur-tp.dt;
            elseif strcmp(trl.expname, 'Troyk_hypothetical')
                trl.order = -1;
                if ~isfield(trl, 'amp');    trl.amp = 1000; end
                if ~isfield(trl, 'pw');     trl.pw = .1*10^-3;     end
                if ~isfield(trl, 'dur');    trl.dur =  50*10^-3;   end% duration in ms
                if ~isfield(trl, 'freq');   trl.freq = 200;        end %NaN if not using a temporal model
                if ~isfield(trl, 'trialdur'); trl.trialdur = 1.7;  end % allows filler at end of trial
                trl.t = 0:tp.dt:trl.dur-tp.dt;
            end
            
            if ~isfield(trl, 'dur');    trl.dur = 1000*10^-3;   end% duration in ms
            if ~isfield(trl, 'trialdur');    trl.trialdur = trl.dur;   end% duration in ms
            trl.t = 0:tp.dt:trl.dur-tp.dt;
            if ~isfield(trl, 'pw');     trl.pw = .25 * 10^-3;      end
            if ~isfield(trl, 'ip');     trl.ip = 0;             end % interphase delay
            if ~isfield(trl, 'lag');    trl.lag = trl.pw;       end% delay before the pulse train begins in ms
            if ~isfield(trl, 'order');  trl.order = 1;          end% 1 = cathodic first, -1  = anodic first
            if ~isfield(trl, 'freq');   trl.freq = 60;          end %NaN if not using a temporal model
            if ~isfield(trl, 'amp');    trl.amp = 100;          end% current amplitude in microAmps
            trl = p2p_c.generate_pt(trl, tp);
            trl.CperTrial = (trl.amp/1000) * trl.dur * trl.freq * trl.pw*10.^3;
            trl.CperPulse = trl.pw * trl.amp/1000;
        end
        
        % generate functions
        function [c, v] = generate_corticalmap(c, v)
            % define cortex meshgrid
            c.x = linspace(.5/c.pixpermm,c.cortexSize(2)-.5/c.pixpermm,c.cortexSize(2)*c.pixpermm) + c.cortexCenter(1) - c.cortexSize(2)/2;
            c.y = linspace(.5/c.pixpermm,c.cortexSize(1)-.5/c.pixpermm,c.cortexSize(1)*c.pixpermm) + c.cortexCenter(2) - c.cortexSize(1)/2;
            [c.X,c.Y] = meshgrid(c.x,c.y);
            
            % Make the orientation and OD maps
            % [pinwheel,OD] = makePinwheelODMaps(x,y,sig)
            sz = size(c.X);
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
            
            % create visual co-ordinates for the entire cortical meshgrid
            c.v2c.z = p2p_c.v2c(c, c.X+sqrt(-1)*c.Y);
            c.v2c.ANG = angle(c.v2c.z) * 180/pi;
            c.v2c.ECC = abs(c.v2c.z);
            c.v2c.X = real(c.v2c.z);
            c.v2c.Y = imag(c.v2c.z);
            c.cropPix = ~(c.v2c.ANG<max(v.angList) & c.v2c.ANG>min(v.angList) & ...
                c.v2c.ECC<=max(v.eccList));
            
            % create a mesh
            v.zAng = linspace(0,max(v.eccList),v.n)'*exp(sqrt(-1)*v.angList*pi/180);
            c.v2c.gridAngZ = p2p_c.c2v(c, v.zAng);
            v.zEcc = [v.eccList'*exp(sqrt(-1)*linspace(-90,90,v.n)*pi/180)]';
            c.v2c.gridEccZ = p2p_c.c2v(c, v.zEcc);
            
            c.RFmap = max(c.slope .* abs(c.v2c.ECC), c.min)+c.intercept;
        end
        
        function [c] = generate_ef(c, varargin)
            % generates an electric field for each electrode, currently assumes they
            % are on the surface
            %
            % max electric field is normalized to 1
            idx = 1:length(c.e);
            
            if ~isfield(c, 'emodel')
                c.emodel = 'Tehovnik';
                I0 = 1;
            end
            for ii=1:length(idx)
                R=sqrt((c.X-c.e(idx(ii)).x).^2+(c.Y-c.e(idx(ii)).y).^2);
                Rd = R-c.e(idx(ii)).radius;
                pt_ef = ones(size(c.X));
                if strcmp(c.emodel, 'Tehovnik')
                    I0=1; k = 6.75; % in cm
                    pt_ef(R>c.e(idx(ii)).radius)=I0./(1+k*Rd(R>c.e(idx(ii)).radius).^2);
                else
                    %    pt_ef(R>c.e(idx(ii)).radius)=2/pi*(asin(c.e(idx(ii)).radius./R(R>c.e(idx(ii)).radius)));
                end
                c.e(idx(ii)).ef = pt_ef;
            end
        end
        function [c,v] = generate_corticalresponse(c, v)
            % generates the sum of weighted receptive fields activated by an electrode
            % normalized so the max is 1
            v.e.rfmap_noRF = zeros(size(v.X)); % percept based on electric field
            v.e.rfmap = zeros([size(v.X), 2]); % percept that includes a cortical model
            
            for pixNum = 1:length(c.X(:))
                if mod(pixNum, 8000)==0
                    disp([num2str(round((100*pixNum)/length(c.X(:)))),  '% complete' ]);
                end
                x0 = c.v2c.X(pixNum); % x center
                y0 = c.v2c.Y(pixNum); % y center
                theta = pi-c.ORmap(pixNum);  %orientation
                sigma_x = c.RFmap(pixNum) * c.ar; % major axis sd
                sigma_y = c.RFmap(pixNum); % minor axis sd
                G = Gauss_2D(v,x0,y0,theta,sigma_x,sigma_y);
                c.target.R(pixNum) = sum(G(:).*v.target.img(:));
            end
            c.target.R = reshape(c.target.R, size(c.X));
            
            function G = Gauss_2D(v,x0,y0,theta,sigma_x,sigma_y)
                % Generates oriented 2D Gaussian on meshgrid v.X,v.Y
                aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
                bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
                cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
                G = exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
            end
        end
        function [v] = generate_rfmap(c, v, varargin)
            % generates the sum of weighted receptive fields activated by an electrode
            % normalized so the max is 1
            idx = 1:length(v.e);
            for ii = 1:length(idx)
                disp([num2str(round((100*ii)/length(idx))),  '% electrodes complete' ]);
                
                
                rfmap_noRF = zeros(size(v.X)); % percept based on electric field
                rfmap = zeros([size(v.X), 2]); % percept that includes a cortical model
                
                for pixNum = 1:length(c.X(:))
                    if c.e(idx(ii)).ef(pixNum) > c.efthr
                        x0 = c.v2c.X(pixNum); % x center
                        y0 = c.v2c.Y(pixNum); % y center
                        if strcmp(c.e(idx(ii)).hemi, 'rh')
                            x0=-x0; y0 = -y0;
                        end
                        theta = pi-c.ORmap(pixNum);  %orientation
                        sigma_x = c.RFmap(pixNum) * c.ar; % major axis sd
                        sigma_y = c.RFmap(pixNum); % minor axis sd
                        
                        % things you need to go from sigmas and theta to G
                        aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
                        bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
                        cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
                        
                        % oriented 2D Gaussian
                        G = c.e(idx(ii)).ef(pixNum) * exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
                        
                        rfmap(:, :, 1)  =   rfmap(:, :, 1)  + c.ODmap(pixNum)*G;
                        rfmap(:, :, 2)  =   rfmap(:, :, 2)  + (1-c.ODmap(pixNum))*G;
                        
                        % scoreboard version
                        rfmap_noRF = rfmap_noRF + (c.e(idx(ii)).ef(pixNum) * exp(-( (v.X-x0).^2/(0.01) + (v.Y-y0).^2/(.01))));
                    end
                end
                if sum(rfmap(:)>0)<20
                    disp('WARNING! Too few pixels passed ef threshold.');
                    disp(' try lowering c.efthr, checking location of electrodes relative to cortical sheet & ');
                    disp('checking the sampling resolution of cortex');
                end
                v.e(idx(ii)).rfmap_noRF = rfmap_noRF./max(rfmap_noRF(:));
                v.e(idx(ii)).rfmap = rfmap./max(rfmap(:));
            end
        end
        function [trl,v] = generate_phosphene(v, tp, trl)
            
            % calculate the neural response over time for each trial
            if ~isnan(trl.freq)
                trl = p2p_c.p2p_finite_element(tp, trl);
                trl.maxphos = v.e(trl.e).rfmap.*max(trl.resp); % the scaling due to current integration over time
            else
                trl.maxphos = v.e(trl.e).rfmap; % the scaling due to current integration
            end
            trl.sim_area = (1/v.pixperdeg.^2) * sum(trl.maxphos(:) > v.drawthr)/2; % mean of left and right eyes
            
            for i=1:2 % left and right eye
                p = p2p_c.fit_ellipse_to_phosphene(trl.maxphos(:,:,i)>v.drawthr,v);
                trl.ellipse(i).x = p.x0;
                trl.ellipse(i).y = p.y0;
                trl.ellipse(i).sigma_x = p.sigma_x;
                trl.ellipse(i).sigma_y = p.sigma_y;
                trl.ellipse(i).theta = p.theta;
            end
            
            % what rule to use to translate phosphene image to brightness?
            
            beta = 6; % soft-max rule across pixels for both eyes
            trl.sim_brightness = ((1/v.pixperdeg.^2) * sum(trl.maxphos(:).^beta)^(1/beta));
        end
        
        function p = fit_ellipse_to_phosphene(img,v)
            M00 = sum(sum(img));
            M10 = sum(sum(v.X.*img));
            M01 = sum(sum(v.Y.*img));
            M11 = sum(sum(v.X.*v.Y.*img));
            M20 = sum(sum(v.X.^2.*img));
            M02 = sum(sum(v.Y.^2.*img));
            
            p.x0 = M10/M00;
            p.y0 = M01/M00;
            
            mu20 = M20/M00 - p.x0^2;
            mu02 = M02/M00 - p.y0^2;
            mu11 = M11/M00 - p.x0*p.y0;
            
            a = (mu20+mu02)/2;
            b = .5*sqrt(4*mu11^2+(mu20-mu02)^2);
            
            lambda_1 = a+b;
            lambda_2 = a-b;
            
            p.theta = -.5*atan2(2*mu11,mu20-mu02);
            p.sigma_x = 2*sqrt(lambda_1);
            p.sigma_y = 2*sqrt(lambda_2);
            
        end
        function  trl = p2p_finite_element(tp, trl )
            % Implements accumulation of current over time using a very simple finite element method
            % written GMB 11/10/2017
            % adapted for cortex IF 3/2/2018
            % Finite difference method:
            % Initial conditions
            
            tmp.chargeacc = 0;
            tmp.ca = 0;
            tmp.cl = 0;
            tmp.R1 = 0;
            tmp.R1 = 0;
            tmp.R4 =  zeros(1,4);
            
            for i=1:length(trl.pt)
                tmp.R1 = tmp.R1 + tp.dt * ((tp.scFac * trl.pt(i))-tmp.R1)/tp.tau1;
                tmp.R2 = max(tmp.R1, 0);
                if ~strcmp(tp.model, 'chronaxie')
                    
                    if strcmp(tp.model, 'nanduri')
                        tmp.R3 = tp.asymptote./(1+exp(-(tmp.R2(:,1)./tp.slope)+(tp.shift)));
                    elseif strcmp(tp.model, 'sigmoid')
                        tmp.R3 = tp.asymptote .* tmp.R2.^2./(tmp.R2.^2 + tp.e50.^2);
                    elseif strcmp(tp.model, 'normcdf')
                        tmp.R3 = tp.asymptote.*normcdf(tmp.R2, tp.mean, tp.sigma);
                    elseif strcmp(tp.model, 'nosaturation')
                        tmp.R3 = tmp.R2;
                    end
                    tmp.R4(:, 1) = tmp.R3;
                    for j=1:3
                        tmp.R4(:, j+1) = tmp.R4(:, j+1) + tp.dt*(tmp.R4(:, j) - tmp.R4(:, j+1))/tp.tau3;
                    end
                    trl.R3(i) = tmp.R3;
                    trl.resp(i)=tmp.R4(:, 4);
                end
                trl.R1(i) = tmp.R1;
                trl.R2(i) = tmp.R2;
            end
        end
        function trl = generate_pt(trl, tp)
            if isnan(trl.freq) % if not using a temporal model at all
                pt = 1;
            else
                t = 0:tp.dt:trl.dur-tp.dt;
                on =  mod(trl.t,1/trl.freq) < trl.pw;
                delay =  trl.pw+trl.ip;
                lag = round(trl.lag/tp.dt);
                off = mod(trl.t-delay,1/trl.freq) < trl.pw;
                tmp  = trl.amp.*(on-off);
                trl.pt= zeros(1, lag+length(tmp));
                trl.pt(lag+1:lag+length(tmp))=tmp;
            end
            if trl.dur<trl.trialdur
                trl.pt((end+1):round((trl.trialdur/tp.dt))) = 0;
                trl.t = 0:tp.dt:1-tp.dt;
            end
        end
        function v = generate_visualtarget(v)
            [x, y] = pol2cart(v.e.ang, v.e.ecc-v.target.offset)
            v.target.img = sqrt(((v.X-x).^2)+((v.Y-y).^2))<v.target.rad;
        end
        
        % transforms
        function c2v_out = c2v(c, z)
            % takes in imaginary numbers, and finds out where cortical values are in visual space (map)
            c2v_out = c.k*log(z + c.a);
        end
        function v2c_out = v2c(c, z)
            % takes in imaginary numbers, places visual values into the cortical grid (mapinv)
            v2c_out = exp(z/c.k)-c.a;
        end
        function logx2raw(base, precision)
            % logx2raw(base, precision)
            %
            % Converts X-axis labels from log to raw values.
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
        
        % plotting functions
        function plotcortgrid(img, c,  varargin)
        % plotcortgrid(img, c)
        % plotcortgrid(img, c, cmap,figNum, evalstr)
        % takes as input:
        %   cortical image
        %   the structure c that defines the cortical surface
        % optional arguments: 
        %   colormap, figure number and a string to evaluate
        %  (e.g. ''title('''corticalsurface''')' or 'subplot(1, 2,1)';
            
            if nargin<3 || isempty(varargin{1});  cmap = gray(256);   else cmap = varargin{1}; end
            if nargin<4 || isempty(varargin{2}); figNum = 1;        else figNum = varargin{2}; end
            if nargin<5 || isempty(varargin{3});  evalstr = '';      else evalstr = varargin{3}; end
            
            
            if isfield(c,'cropPix')
                img(c.cropPix) = NaN;
                img= img+2;
                cmap = [0,0,0;cmap];
            end
            
            fH=figure(figNum);
            eval(evalstr); colormap(cmap);
            if ~isempty(img)
                image(c.x, c.y, img); hold on
            end
            xlabel('mm'); ylabel('mm')
            set(gca,'YDir','normal');
            plot(c.v2c.gridAngZ, '-', 'Color', c.gridColor);
            plot(c.v2c.gridEccZ, '-', 'Color', c.gridColor);
            
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

            if nargin<3 || isempty(varargin{1});  cmap = gray(256);   else cmap = varargin{1}; end
            if nargin<4 || isempty(varargin{2});  figNum = 1;        else figNum = varargin{2}; end
            if nargin<5 || isempty(varargin{3});  evalstr = '';      else evalstr = varargin{3}; end
            
        
            fH=figure(figNum); hold on
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
        function figure_xyplot(trl, xval, yval1, yval2, figNum, substr, varargin)
            figure(figNum); hold on
            eval(substr)
            
            if nargin <7
                opts.marker = 'o';
                opts.color = 'k';
                opts.MarkerSize = 8;
                opts.LineStyle = 'none';
            else
                opts = varargin{1};
            end
            
            if ~isfield(opts, 'marker'); opts.marker = 'o'; end
            if ~isfield(opts, 'color'); opts.color = [0 0 0]; end
            if ~isfield(opts, 'MarkerSize'); opts.MarkerSize = 8; end
            if ~isfield(opts, 'LineStyle'); opts.LineStyle = 'none'; end
            set(gca, 'FontSize', 8); hold on
            
            x  = eval(['cat(1,trl(:).', xval,');']);
            y1 = eval(['cat(1,trl(:).', yval1,');']);
            if ~isempty(yval2)
                y2 = eval(['cat(1,trl(:).', yval2,');']);
            end
            switch xval
                case 'CperTrial'
                    x = log(x);
                    logFlag_x = 1;
                    xlabel('Charge Deposited per Trial (?C)');
                case 'CperPulse'
                    logFlag_x = 0;
                    xlabel('Charge Deposited per Pulse (?C)');
                case 'amp'
                    logFlag_x = 0;
                    xlabel('Current');
                case 'draw_brightness'
                    logFlag_x = 0;
                    xlabel('Phosphene brightness');
                case 'draw_area'
                    x = log(x);
                    logFlag_x = 1;
                    xlabel('Phosphene area (deg^2)');
                case  'draw_radius'
                    logFlag_x = 0;
                    xlabel('Phosphene radius (deg)');
                case 'draw_diameter'
                    logFlag_x = 0;
                    xlabel('Phosphene diameter (deg)');
                otherwise
                    logFlag_x = 0;
                    xlabel('xval')
            end
            
            switch yval1
                case 'draw_brightness'
                    y2=y2*[nanmean(y1)./nanmean(y2)];
                    logFlag_y = 0;
                    ylabel('Phosphene brightness')
                case 'draw_area'
                    y1 = log(y1);% y2 = log(y2);
                    logFlag_y = 1;
                    ylabel('Phosphene area (deg^2)')
                case 'draw_radius'
                    logFlag_y = 0;
                    ylabel('Phosphene radius (deg)')
                case 'sim_area'
                    y1 = log(y1);
                    logFlag_y = 1;
                    ylabel('SIM Phosphene area (deg^2)')
                case 'sim_radius'
                    logFlag_y = 0;
                    ylabel('SIM Phosphene radius (deg)')
                case 'sim_diameter'
                    logFlag_y = 0;
                    ylabel('Phosphene diameter (deg)')
                case 'norm_sim_radius'
                    y1 = y1./max(y1(:));
                    logFlag_y = 0;
                    ylabel('Phosphene radius (deg)')
                case 'sim_brightness'
                    y1=y1*[nanmean(x)./nanmean(y1)];
                    logFlag_y = 0;
                    ylabel('SIM Phosphene brightness')
                otherwise
                    logFlag_ = 0;
                    ylabel('yval')
            end
            [~, idx]=sort(x);
            idx = idx(y1(idx)>0);
            
            plot(x(idx), y1(idx), 'Marker', opts.marker, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', opts.color, 'MarkerSize', opts.MarkerSize, 'LineStyle', opts.LineStyle); hold on
            if ~isempty(yval2)
                plot(x(idx), y2(idx), 'Marker', opts.marker, 'MarkerEdgeColor', opts.color, 'MarkerFaceColor', 'none', 'MarkerSize', opts.MarkerSize, 'LineStyle', opts.LineStyle);
            end
            
            if logFlag_x
                set(gca, 'XLim', log([1 2000]));
                set(gca, 'XTick', log([1 10 100 1000]));
                logx2raw;
            end
            if logFlag_y
                set(gca, 'YLim', log([.005 500]));
                set(gca, 'YTick', log([.1 1 100]));
                logy2raw;
            end
            if isempty(yval2)
                switch xval
                    case 'draw_area'
                        set(gca, 'XLim', log([.005 500]));
                        set(gca, 'XTick', log([.1 1 100]));
                    case 'draw_radius'
                        axis([0 35 0 35]);
                    case 'draw_brightness'
                        axis([0 25 0 25]);
                end
            end
        end
        
        function p2p_generic()
            c = define_cortex();
            v = define_visualmap();
            c = generate_corticalmap(c, v);
            
            c = define_electrodes(c, v);
            
            tp = define_temporalparameters();
            trls = define_exp(tp);
            
            plotcortgrid(64 * (c.ORmap+pi)/(pi*2), c, 'Orientation pinwheels', hsv(64), 1, 'subplot(1, 1, 1)');
            plotcortgrid(64 * c.ODmap, c, 'Ocular dominance columns', gray(64), 2, 'subplot(1, 1, 1)');
            plotcortgrid(15* c.RFmap, c, 'Receptive field size', hot(64), 3, 'subplot(1, 1, 1)');
            
            c = generate_ef(c);
            plotcortgrid(c.e.ef * 64, c, 'electric field', gray(64), 4, 'subplot(1,1,1)');
            v = generate_rfmap(c, v);
            plotretgrid(v.e.rfmap(:, :, 1)*64, c, v,'rfmap', gray(64), 10, 'subplot(1,3,1)');
            plotretgrid(v.e.rfmap_noRF*64, c, v, 'rfmap no RF', gray(64), 10, 'subplot(1,3,2)');
            
            v = generate_phosphene(v, tp, trls);
            plotretgrid(v.trls(1).maxphos(:, :, 1)*64, c, v,'phosphene', gray(64), 10, 'subplot(1,3,3)');
        end
        
        function junk()
            %% CODE TO INTEGRATE
            %         idx=find(trls(t).pt>0);
            %         if ii==5; idx=setdiff(idx,3:5); end % drawings not obtained
            %
            %         poly_area = D(ii).poly_area.val(idx);
            %         poly_area(isnan(poly_area))=0;
            %
            %         if ~c.smallflag
            %             p.thr = 50;  p=UWfit('findDrawingThreshold', p, {'thr'}, Dsim(ii).EF.normalizedphosphene(:), Dsim(ii).CI.val(idx), ...
            %                 (poly_area * c.pixperdeg^2));
            %             Dsim(ii).EF.drawingthreshold = p.thr;
            %             disp(['EF drawing threshold = ', num2str(Dsim(ii).EF.drawingthreshold)]);
            %
            %             p.thr = 50; p=UWfit('findDrawingThreshold', p, {'thr'}, Dsim(ii).CM.normalizedphosphene(:), Dsim(ii).CI.val(idx), ...
            %                 (poly_area * c.pixperdeg^2));
            %             Dsim(ii).CM.drawingthreshold = p.thr;
            %             disp(['CM drawing threshold = ', num2str(Dsim(ii).CM.drawingthreshold)]);
            %         end
            %         for cc=1:length(idx)
            %             Dsim(ii).EF.poly_area.val(idx(cc)) = (1/(c.pixperdeg^2)) * ...
            %                 length( find((Dsim(ii).EF.normalizedphosphene(:) * Dsim(ii).CI.val(idx(cc)) ) ...
            %                 >  Dsim(ii).EF.drawingthreshold) );
            %             Dsim(ii).CM.poly_area.val(idx(cc)) = (1/(c.pixperdeg^2)) * ...
            %                 length( find((Dsim(ii).CM.normalizedphosphene(:) *  Dsim(ii).CI.val(idx(cc)) ) ...
            %                 >  Dsim(ii).CM.drawingthreshold) );
            %         end
            %
            % end
            
            %% NOTES
            % USING CORTICAL POSITION
            % Winawer site one is LH, all the rest are RH
            % original values based on the Benson template extracted from Winawer, are converted into x y co-ordinates using:
            % Z=eData(s,1).*exp(sqrt(-1)*eData(s,2)*pi/180); xe=real(Z); ye=imag(Z);
            % orig values in visual space co-ordinates were ecc, angle =
            % ea_ret = [7.8774	8.8141; 18.2743	23.5752; 2.4848	48.3118 ; 7.8671	166.642 ;  3.337	49.077]; % pulled from Win paper, in degrees
            % % convert so all projecting to the same hemisphere
            % ea_ret = [7.8774	8.8141; 18.2743	23.5752; 2.4848	48.3118 ; 7.8671	180-166.642 ;  3.337	49.077]; % ecc, angle, pulled from Win paper, in degrees
            
            % USING PHOSPHENE LOCATION
            %xyloc = [25 9; -7.65 -1.85; -4.05 3.14; -1.35 1.35; -.75 .5]; % pulled using ginput from fig 2
            % ea_ret = [ 26.6 19.8; 9  -166.4; 5.12 142.2; 1.9 135; 0.90 146.3]; % ecc, angle
            % ea_ret = [ 26.6 19.8; 9  -166.4; 5.12 142.2; 1.9 135; 1 146.3]; % ecc, angle
            
            %ea_ret = [ 26.6 19.8; 9  -(180-166.4); 5.12 (180-142.2); 1.9 (180-135); 1 (180-146.3)]; % ecc, angle, flipped 2:5 to right visual
        end
    end
end

