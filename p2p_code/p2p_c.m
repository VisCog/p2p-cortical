% function [c,trl] = p2p_main()


classdef p2p_c
    methods(Static)
        
        %p2p_Winawer();
        
        %p2p_Tehovnik();
        %p2p_generic();
        
        
     
        % definitions
        function v = define_visualmap(v)
            if ~isfield(v.e, 'ang')
                v.e.ang = 19.8; % position of the electode in visual co-ordinates
                v.e.ecc = 26.6;
            end
            if ~isfield(v, 'retinaSize')
                v.retinaSize = [70,70]; %  [height, width diameter in degrees]
            end
            
            if ~isfield(v,'retinaCenter')
                v.retinaCenter = [0,0];
            end
            
            if ~isfield(v,'pixperdeg')
                v.pixperdeg = 20;
            end
            if ~isfield(v, 'drawthr')
                v.drawthr = 0.15;
            end
            v.x = linspace(.5/v.pixperdeg,v.retinaSize(2)-.5/v.pixperdeg,v.retinaSize(2)*v.pixperdeg)+v.retinaCenter(1) - v.retinaSize(2)/2;
            v.y = linspace(.5/v.pixperdeg,v.retinaSize(1)-.5/v.pixperdeg,v.retinaSize(1)*v.pixperdeg)+v.retinaCenter(2) - v.retinaSize(1)/2;
            [v.X,v.Y] = meshgrid(v.x, v.y);
            
            %Make the grid in retinal coordinates
            v.angList = -90:45:90;
            v.eccList = [1 2 3 5 8 13 21 34];
            v.gridColor = [1 1 0];
            v.n = 201;
        end
        function c = define_cortex(c)
            
            % cortical magnification, typical log z transformation parameters (Based on Duncan and Boynton)
            c.k = 20; %scale
            if ~isfield(c, 'a')
                c.a = 0.5; %fovea expansion for human, macaque is 0.3
            end
            
            % current spread parameters
            % current spread parameters, based on Ahuja
            %             c.afac=1.69; % parameters for current spread
            %             c.ct=14;
            
            % receptive field parameters
            c.ar = .5;  % aspect ratio of V1 RFs
            c.sig = .5; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex. sig determines the distribution of OD values. Default is 5.  The larger sig, the more the distribution tends toward 0 and 1.
            c.ODsize = 0.863; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
            c.filtSz = 3; % 3mm creates the initial OD and orientation maps
            
            if ~isfield(c, 'slope')
                c.slope = .1652;
                c.min = 0.1;
                c.intercept = .1;
            end
            
            % define the size and resolution of cortical and visual space parameters
            c.gridColor = [1,1,0];
            if ~isfield(c, 'cortexSize')
                c.cortexSize = [80,100];  % %[height, width] Size of cortical maps (mm)
            end
            if ~isfield(c, 'cortexCenter')
                c.cortexCenter = [30,0]; % center of electrode array (mm on cortex)
            end
            if ~isfield(c, 'pixpermm')
                c.pixpermm = 8; % choose the resolution to sample in mm.
            end
        end
        function c = define_electrodes(c, v, varargin)
            % either needs an electrode position in cortical co-ordinates
            % or needs to take in the position of the electrode in visual co-ordinates
            if nargin <3
                idx = 1:length(v.e);
            else
                idx = varargin{1};
            end
            
            if ~isfield(c.e, 'radius')
                for ii=1:length(idx)
                    c.e(idx(ii)).radius = 500/1000;
                end
            end
            
            for ii=1:length(idx)
                if v.e(idx(ii)).ang>-90 && v.e(idx(ii)).ang <90
                    c.e(idx(ii)).hemi = 'lh';
                else
                    c.e(idx(ii)).hemi = 'rh';
                end
                c.e(idx(ii)).area = pi*(c.e(idx(ii)).radius.^2);
            end
            
            % tranform electrodes
            if isfield(v.e, 'ecc')
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
            tp.dt = .01 * 10^-3; % time sampling in ms
            if ~isfield(tp, 'scFac'); tp.scFac = 1; end
            tp.tau1 = .2 * 10^-3; %Tehovnik et al 2004
            tp.tau2_ca = 45.250* 10^-3;  %38-57, from retina
            tp.tau3 =  26.250* 10^-3; % 24-33 from retina
            tp.e = 8.73; % from retina
            
            % leak out of charge accumulation
            tp.flag_cl=0; % 1 if you want to charge to leak back out of the system
            tp.tau2_cl = tp.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca
            
            % nonlinearity parameters
            tp.slope=.3; % larger number = shallower slope
            tp.asymptote=14;
            tp.shift=47; % shifts curve along x-axis
        end
        function trl = define_trial(tp,trl)
            if ~isfield(trl,'expname'); trl.expname = 'generic'; end
            if ~isfield(trl,'e'); trl.e = 1; end
            
            if strcmp(trl.expname, 'Tehovnik')
                trl.dur = 100*10^-3; % duration in ms
                trl.pw = .2/1000; % pulse width in ms
                trl.order = -1; % 1 = cathodic first, -1  = anodic first
                trl.freq = 200; % -1 for a single pulse, NaN if not using a temporal model
                trl.amp = 50; % current amplitude in microAmps
            elseif strcmp(trl.expname, 'Bosking')
                trl.dur = 250*10^-3; % duration in ms
                trl.t = 0:tp.dt:trl.dur-tp.dt;
                trl.pw = .1/1000; % pulse width in ms
                trl.order = 1; % 1 = cathodic first, -1  = anodic first
                trl.freq = 200; % -1 for a single pulse, NaN if not using a temporal model
            elseif strcmp(trl.expname,'Evans')
                trl.order = -1;
            if ~isfield(trl, 'pw');     trl.pw = .25/1000;     end
            if ~isfield(trl, 'dur');    trl.dur = .5;   end% duration in ms
            if ~isfield(trl, 'freq');   trl.freq = 50;          end %NaN if not using a temporal model
                trl.t = 0:tp.dt:trl.dur-tp.dt;
            end

            if ~isfield(trl, 'dur');    trl.dur = 1000*10^-3;   end% duration in ms
            trl.t = 0:tp.dt:trl.dur-tp.dt;
            if ~isfield(trl, 'pw');     trl.pw = .25 * 10^-3;      end
            if ~isfield(trl, 'ip');     trl.ip = 0;             end % interphase delay
            if ~isfield(trl, 'lag');    trl.lag = trl.pw;       end% delay before the pulse train begins in ms
            if ~isfield(trl, 'order');  trl.order = 1;          end% 1 = cathodic first, -1  = anodic first
            if ~isfield(trl, 'freq');   trl.freq = 60;          end %NaN if not using a temporal model
            if ~isfield(trl, 'amp');    trl.amp = 100;          end% current amplitude in microAmps
            trl.pt = p2p_c.generate_pt(trl, tp);
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
            if nargin <2
                idx = 1:length(c.e);
            else
                idx = varargin{1};
            end
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
        function [c] = generate_corticalresponse(c, v)
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
            
            if nargin <3
                idx = 1:length(v.e);
            else
                idx = varargin{1};
            end
            
            for ii = 1:length(idx)
                rfmap_noRF = zeros(size(v.X)); % percept based on electric field
                rfmap = zeros([size(v.X), 2]); % percept that includes a cortical model
                
                for pixNum = 1:length(c.X(:))
                    if mod(pixNum, 12000)==0
                        disp([num2str(round((100*pixNum)/length(c.X(:)))),  '% complete' ]);
                    end
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
                if sum(rfmap(:)>0)<100
                    disp('WARNING! Too few pixels passed ef threshold, try lowering c.efthr');
                end
                v.e(idx(ii)).rfmap_noRF = rfmap_noRF./max(rfmap_noRF(:));
                v.e(idx(ii)).rfmap = rfmap./max(rfmap(:));
            end
        end
        function trl = generate_phosphene(v, tp, trl)
            
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
            
            tmp.R3norm = 0;
            tmp.R4a =  zeros(1,4);
            %tmp.R4ax =  zeros(1,4);

            for i=1:length(trl.pt)-1
                % R1
                tmp.R1= tmp.R1 + tp.dt * ((tp.scFac * trl.pt(i))-tmp.R1)/tp.tau1;
                tmp.R4a(:,1) = max(tmp.R1, 0);
                tmp.R4a(:, 1)=tp.asymptote./(1+exp(-(tmp.R4a(:,1)./tp.slope)+tp.shift));
                for j=1:3
                    tmp.R4a(:,j+1) = tmp.R4a(:,j+1) + tp.dt*(tmp.R4a(:,j) - tmp.R4a(:,j+1))/tp.tau3;
                end
                tmp.R4(i) = tmp.R4a(:,4);
            end
            trl.resp=tmp.R4;
        end
        function pt = generate_pt(trl, tp)
            if isnan(trl.freq) % if not using a temporal model at all
                pt = 1;
            else
                t = 0:tp.dt:trl.dur-tp.dt;
                on =  mod(t,1/trl.freq) < trl.pw;
                delay =  trl.pw+trl.ip;
                lag = round(trl.lag/tp.dt);
                off = mod(t-delay,1/trl.freq) < trl.pw;
                tmp  = trl.amp.*(on-off);
                pt= zeros(1, lag+length(tmp));
                pt(lag+1:lag+length(tmp))=tmp;
                pt=pt(1:length(tmp));
            end
            if trl.dur<1
                pt((end+1):round((1/tp.dt))) = 0;
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
        function plotcortgrid(img, c,  cmap,figNum, spstr)
            
            if isfield(c,'cropPix')
                img(c.cropPix) = NaN;
                img= img+2;
                cmap = [0,0,0;cmap];
            end
            
            fH=figure(figNum);
            eval(spstr); colormap(cmap);
            image(c.x, c.y, img); hold on
            xlabel('mm'); ylabel('mm')
            set(gca,'YDir','normal');
            plot(c.v2c.gridAngZ, '-', 'Color', c.gridColor);
            plot(c.v2c.gridEccZ, '-', 'Color', c.gridColor);
            
            axis equal;  axis tight
            set(gca,'XLim',[min(c.x(:)),max(c.x(:))]);
            set(gca,'YLim',[min(c.y(:)),max(c.y(:))]);
            drawnow;
        end
        function plotretgrid(img, v, cmap, figNum, spstr)
            
            fH=figure(figNum); hold on
            eval(spstr);
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
                    xlabel('Charge Deposited per Trial (µC)');
                case 'CperPulse'
                    logFlag_x = 0;
                    xlabel('Charge Deposited per Pulse (µC)');
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
        
        function p2p_Tehovnik()
            c.a = 0.3; c.cortexSize = [20,45]; c.pixpermm = 10; c.efthr = .15;
            c.slope = .051; c.min = 0; c.intercept = .214;
            c = define_cortex(c);
            v.e.ang = 0; v.e.ecc = 3.7;v.retinaSize = [5,15];v.pixperdeg = 40; v = define_visualmap(v);
            c = generate_corticalmap(c, v);
            c.e.radius = 50/1000; c = define_electrodes(c, v);
            
            v.target.offset = 0.5; v.target.rad = 0.2; v = generate_visualtarget(v);
            c = generate_corticalresponse(c, v);
            c = generate_ef(c);
            v = generate_rfmap(c, v);
            tp = define_temporalparameters();
            trl = define_exp(tp, 'Tehovnik');
            v = generate_phosphene(v, tp, trls);
            
            % cortical projection into visual field
            plotretgrid(v.e.rfmap(:, :, 1)*1000 + v.target.img*64, c, v,'rfmap', gray(64), 10, 'subplot(1,3,1)');
            plotretgrid(v.e.rfmap_noRF*64 + v.target.img*64, c, v, 'rfmap no RF', gray(64), 10, 'subplot(1,3,2)');
            plotretgrid(v.trls(1).maxphos(:, :, 1)*1000 + v.target.img*64, c, v,'phosphene', gray(64), 10, 'subplot(1,3,3)');
            
            % plotcortgrid(c.e.ef * 64, c, 'electric field', gray(64), 4, 'subplot(2,2,1)');
            plotcortgrid(4 * c.target.R, c, 'Cortical response to visual target', gray(64), 4, 'subplot(1, 2, 1)');
            plotcortgrid(64 * c.e.ef, c, 'Electrical response', gray(64), 4, 'subplot(1, 2, 2)');
            
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

    