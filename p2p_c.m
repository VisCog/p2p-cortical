% p2p_c
%
%  holds all support functions for p2p cortex
% project.
%
% functions can be called from outside with 'p2p_c.<function name>'


classdef p2p_c
    methods(Static)
        %% definitions - temporal properties

  
        function trl = define_trial(tp, varargin)
            if nargin < 2;  trl = [];
            else; trl = varargin{1}; end

            if ~isfield(trl,'e'); trl.e = 1; end

            % durations are all in s
            if ~isfield(trl, 'dur');    trl.dur = 1000*10^-3;   end % duration in s of the electrical stimulation
            if ~isfield(trl, 'simdur');    trl.simdur = trl.dur+250*10^-3;   end
            % duration in s of the length of time being simulated, needs to be
            % longer to allow time for the neural response

            trl.t = 0:tp.dt:trl.dur-tp.dt;
            if ~isfield(trl, 'pw');     trl.pw = .1 * 10^-3;    end
            if ~isfield(trl, 'ip');     trl.ip = 0;             end % interphase delay
            if ~isfield(trl, 'lag');    trl.lag = 2*trl.pw;       end % delay before the pulse train begins
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
                    off = mod(trl.t-trl.pw,1/trl.freq) <=trl.pw;
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
      

        %% cortical stuff
        function c = define_cortex(c)
            % cortical magnification,
            % typical log z transformation parameters (based on early
            % Schwartz model
            if ~isfield(c, 'animal');   c.animal = 'human'; end
            if ~isfield(c, 'efthr'); c.efthr = 0.05; end % what mag of electric field goes through the model, just a  speed thing
            if strcmp(c.animal, 'human')
                c.k = 15; %scale
                c.a = 0.5; %fovea expansion for human, macaque is 0.3
            elseif strcmp(c.animal, 'macaque')
                c.k = 5; %scale
                c.a = 0.3; % values set by eyeballing Toottell data
            end
            % receptive fields
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
            % ocular dominance columns
            if ~isfield(c, 'sig'); c.sig = .5;  end
            % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex. sig determines the distribution of OD values. Default is 5.
            % The larger sig, the more the distribution tends toward 0 and 1.

            if strcmp(c.animal, 'human')
                c.ODsize = 0.863; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 3; % 3mm creates the initial OD and orientation maps
            else
                c.ODsize = 0.531; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 1.85; % 3mm creates the initial OD and orientation maps
            end

            % define the size and resolution of cortical and visual space parameters
            c.gridColor = [1,1,0];
            if ~isfield(c, 'cortexSize'); c.cortexSize = [80,100]; end %[height, width] Size of cortical maps (mm)
            if ~isfield(c, 'cortexCenter'); c.cortexCenter = [30,0]; end% center of electrode array (mm on cortex)
            if ~isfield(c, 'pixpermm');  c.pixpermm = 8; end    % choose the resolution to sample in mm.
        end
        function [c] = generate_ef(c, varargin)
            % generates an electric field for each electrode, currently assumes they
            % are on the surface

            % max electric field is normalized to 1
            idx = 1:length(c.e);

            if ~isfield(c, 'emodel')
                c.emodel = 'Tehovnik';      I0 = 1;   end
            for ii=1:length(idx)
                R=sqrt((c.X-c.e(idx(ii)).x).^2+(c.Y-c.e(idx(ii)).y).^2);
                Rd = R-c.e(idx(ii)).radius;
                pt_ef = ones(size(c.X));
                if strcmp(c.emodel, 'Tehovnik')
                    I0=1; k = 6.75; % in cm
                    pt_ef(R>c.e(idx(ii)).radius)=I0./(1+k*Rd(R>c.e(idx(ii)).radius).^2);
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

            c.RFmap = max(c.slope .* abs(c.v2c.ECC), c.min) + c.intercept;
        end

        function v = define_visualmap(v)
            if ~isfield(v, 'retinaSize');    v.retinaSize = [70,70]; end%  [height, width diameter in degrees]
            if ~isfield(v,'retinaCenter');    v.retinaCenter = [0,0];  end

            if ~isfield(v,'pixperdeg');     v.pixperdeg = 7;       end
            if ~isfield(v, 'drawthr');     v.drawthr = 0.55;    end
            v.x = linspace(.5/v.pixperdeg,v.retinaSize(2)-.5/v.pixperdeg,v.retinaSize(2)*v.pixperdeg)+v.retinaCenter(1) - v.retinaSize(2)/2;
            v.y = linspace(.5/v.pixperdeg,v.retinaSize(1)-.5/v.pixperdeg,v.retinaSize(1)*v.pixperdeg)+v.retinaCenter(2) - v.retinaSize(1)/2;
            [v.X,v.Y] = meshgrid(v.x, v.y);

            %Make the grid in retinal coordinates
            if ~isfield(v, 'angList');   v.angList = -90:45:90;    end
            if ~isfield(v, 'eccList');  v.eccList = [1 2 3 5 8 13 21 34];    end

            v.gridColor = [1 1 0];
            v.n = 201;
        end
        function c = define_electrodes(c, v)
            %  takes in the position of the electrode in visual
            %  co-ordinates and pop them onto the cortical surface
            idx = 1:length(v.e);
            if ~isfield(c.e, 'radius')
                for ii=1:length(idx);        c.e(idx(ii)).radius = 500/1000;       end
            end
            if ~isfield(c.e, 'shape')
                for ii=1:length(idx);     c.e(idx(ii)).shape = 'round';    end
            end

            for ii = 1:length(idx)
                c.e(idx(ii)).area = pi*(c.e(idx(ii)).radius.^2);
            end
            for ii=1:length(idx)
                if v.e(idx(ii)).ang>-90 && v.e(idx(ii)).ang <90
                    c.e(idx(ii)).hemi = 'rh';
                else
                    c.e(idx(ii)).hemi = 'lh';
                end
            end

            % transform electrodes
            for ii=1:length(idx)
                if strcmp(c.e(idx(ii)).hemi, 'lh')
                    z = v.e(idx(ii)).ecc.*exp(sqrt(-1)*(v.e(idx(ii)).ang+180)*pi/180);
                else
                    z = v.e(idx(ii)).ecc.*exp(sqrt(-1)*v.e(idx(ii)).ang*pi/180);
                end
                c.e(idx(ii)).z = p2p_c.c2v(c,z);
                c.e(idx(ii)).x = real(c.e(idx(ii)).z);
                c.e(idx(ii)).y = imag(c.e(idx(ii)).z); % turn into mm
            end

        end
        function [v] = generate_rfmap(c, v)
            % calculates spatial phosphenes in visual space based on electrodes in
            % cortical space
            if ~isfield(c, 'rftype')
                c.rftype = 'rf';
            end
            % generates the sum of weighted receptive fields activated by an electrode
            % normalized so the max is 1
            idx = 1:length(v.e);
            for ii = 1:length(idx)
                disp([num2str(round((100*ii)/length(idx))),  '% electrodes complete' ]);
                rfmap = zeros([size(v.X), 2]); % percept that includes a cortical model

                for pixNum = 1:length(c.X(:))
                    if c.e(idx(ii)).ef(pixNum) > c.efthr * 255
                        x0 = c.v2c.X(pixNum); % x center
                        y0 = c.v2c.Y(pixNum); % y center
                        if strcmp(c.e(idx(ii)).hemi, 'lh')
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
                        if strcmp(c.rftype, 'scoreboard')
                            % scoreboard version
                            tmp = double( c.e(idx(ii)).ef(pixNum))/255;
                            rfmap(:, :, 1) = rfmap(:, :, 1) + (tmp * exp(-( (v.X-x0).^2/(0.01) + (v.Y-y0).^2/(.01))));
                            rfmap(:, :, 2) = rfmap(:, :, 1);
                        else
                            tmp = double(c.e(idx(ii)).ef(pixNum))/255;
                            G = tmp * exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
                            rfmap(:, :, 1)  =   rfmap(:, :, 1)  + c.ODmap(pixNum)*G;
                            rfmap(:, :, 2)  =   rfmap(:, :, 2)  + (1-c.ODmap(pixNum))*G;
                        end
                    end
                end
                if sum(rfmap(:)>0)<20
                    disp('WARNING! Too few pixels passed ef threshold.');
                    disp(' try lowering c.efthr, checking location of electrodes relative to cortical sheet & ');
                    disp('checking the sampling resolution of cortex');
                end

                v.e(idx(ii)).rfmap = rfmap./max(rfmap(:));
            end
        end
        function [trl,v] = generate_phosphene(v, tp, trl)
            % calculate the neural response over time across the spatial
            % maps generated by generate_rfmap

            if isnan(trl.freq); trl.maxphos = v.e(trl.e).rfmap; trl.resp = 1;% the scaling due to current integration
            else
                if ~isfield(trl, 'resp') % calculate the response to the pulse train if not already calculated
                    trl = p2p_c.convolve_model(tp, trl);
                end
                trl.maxphos = v.e(trl.e).rfmap.*max(trl.resp); %  scaling the phosphene based on neural response over time
            end
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
            if ~isfield(tp, 'tau1');   tp.tau1 =0.0003; % fixed based on Nowak and Bullier, 1998
            end
            if ~isfield(tp, 'tSamp');   tp.tSamp = 1000;end % subsampling to speed things up
            % fit based on Brindley, 1968, Tehovnk 2004 estimates 0.13-0.24 ms
            if ~isfield(tp, 'tau2');   tp.tau2 =  0.2568;  end     % 24-33 from retina
            if ~isfield(tp, 'ncascades');  tp.ncascades = 3;   end% number of cascades in the slow filter of the temporal convolution
            if ~isfield(tp, 'gammaflag');   tp.gammaflag = 1;   end            %  include second stage game
            

            % leak out of charge accumulation
            %   tp.flag_cl=0; % 1 if you want to charge to leak back out of the system
            %   tp.tau2_cl = tp.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca

            % nonlinearity parameters
            if ~isfield(tp, 'model');  tp.model = 'compression'; end
            if ~isfield(tp, 'power'); tp.power = 20.15;
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
                clear trl;  trl.pw = T.pw(i);   trl.amp = 1;    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 1; %sec
                trl = p2p_c.define_trial(tp,trl);

                tp.scFacInd = 1;     % separate 'scFac' for each experiment
                if isfield(tp,'experimentList')  % set tp_thresh accordingly for this trial
                    experimentNum = find(strcmp(T.experiment{i},tp.experimentList));
                    if ~isempty(experimentNum)  % set thresh_resp for this experiment.
                        tp.scFacInd = experimentNum;
                    end
                end
                thresh(i)= p2p_c.find_threshold(trl,tp);
            end

            err = nansum((thresh-T.amp').^2);
            %    disp(sprintf('tau1 = %g, tau2 = %g, power = %5.2f err= %5.4f  scFac= %5.4f',  tp.tau1,tp.tau2,tp.power,err, tp.scFac));
            disp(fprintf('mean = %g, sigma = %g,  err= %5.4f\n',  tp.mean,tp.sigma,err));
            if isfield(tp,'experimentList')
                for i = 1:length(tp.experimentList)
                    disp(fprintf('%10s: %g\n',tp.experimentList{i},tp.scFac(i)));
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
                clear trl;  trl.pw = T.pw(i);   trl.amp = 1;    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 1; %sec
                trl = p2p_c.define_trial(tp,trl);
                thresh(i)= p2p_c.find_threshold(trl,tp);
            end
            
            err = nansum((thresh-T.amp').^2);
            disp(sprintf('err= %5.4f',  err));

            if isfield(tp,'experimentList')
                for i = 1:length(tp.experimentList)
                    disp(sprintf('%10s: %g',tp.experimentList{i},tp.scFac(i)));
                end
            end
        end
        function err = fit_brightness(tp, T)
            [loop_trl] = p2p_c.loop_convolve_model(tp,T);
            y_est = [loop_trl.maxresp]; y = [T.brightness];
            y_est = reshape(y_est, length(y), 1);
            y = reshape(y, length(y), 1);
            ind = ~isnan(y_est) & ~isnan(y);
            err = -corr(y(ind), y_est(ind));
            disp(sprintf('tau1 =%5.4f, tau2 =%5.4f, power =%5.4f, corr = %5.4f',  tp.tau1, tp.tau2, tp.power,-err));
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
                clear trl;  trl.pw = T.pw(i);   trl.amp = T.amp(i);    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 1; %sec
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
            spikeId = zeros(size(Rtmp));
            for i=1:(length(ptid)-1)
                tNow = trl.t(ptid(i+1));
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
                    lastSpikeTime = trl.t(ptid(i+1));
                else
                    wasRising =1;
                end
            end
            R1= Rtmp*1000;
            R1(R1<0)=0;
            trl.spikes = R1;
            if tp.gammaflag
                if ~isfield(trl, 'imp_resp')
                    dt = t(2)-t(1);
                    h = p2p_c.gamma(tp.ncascades,tp.tau2,t);            % Generate the n-cascade impulse response
                    tid = find(cumsum(h)*dt>.999,1,'first'); % Shorten the filter if needed to speed up the code.
                    if ~isempty(tid)
                        h = h(1:tid);
                    else
                        sprintf('Warning: gamma hdr might not have a long enough time vector');
                    end
                    trl.imp_resp = h;  % close enough to use h
                end
                impFrames = [0:(length(trl.imp_resp)-1)];
                resp = zeros(1,length(t)+length(trl.imp_resp));        % zero stuff out
                for i=1:length(trl.spikes)
                    if spikeId(i)
                        id = find(t>trl.t(ptid(i)),1,'first');
                        resp(id+impFrames)  =   ...
                            resp(id+impFrames) + trl.spikes(i)*trl.imp_resp;
                    end
                end
                resp = p2p_c.nonlinearity(tp, resp);
                trl.maxresp = max(resp); % detection when maxresp goes above a threshold
            else
                resp= trl.spikes;
                trl.maxresp = max(resp);
            end
            % save the time-course of the response for output
            trl.resp = resp(1:length(t));
            trl.tt = t;  % for plotting
        end


        %% utilities
        function out=chronaxie(p,pw)
            out = p.amp./(p.tau*(1-exp(-pw/p.tau)));
        end

        function y=gamma(n,k,t)
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
            if ~isfield(tp, 'scFac');    scFac = 1;
            else scFac  =  tp.scFac; end
            % some of our favorite static nonlinearities:
            switch tp.model
                case 'sigmoid'
                    y = scFac .* x.^tp.power./(x.^tp.power + tp.sigma.^2);
                case 'normcdf'
                    y = normcdf(x, tp.mean, tp.sigma);
                    y(y<0) = 0;
                case 'weibull'
                    y = scFac*p2p_c.weibull(tp,x);
                case 'power'
                    y = scFac*x.^tp.power;
                case 'exp'
                    y = scFac*x.^tp.k;
                case 'compression'
                    y = tp.power.*tanh(x/tp.power);
                case 'linear'
                    y = x;
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

        %% transforms
        % transforms
        function c2v_out = c2v(c, z)
            % takes in imaginary numbers, and finds out where cortical values are in visual space (map)
            c2v_out = c.k*log(z + c.a);
        end
        function v2c_out = v2c(c, z)
            % takes in imaginary numbers, places visual values into the cortical grid (mapinv)
            v2c_out = exp(z/c.k)-c.a;
        end

        %% plotting functions
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

            if nargin<3 || isempty(varargin{1});  cmap = gray(256);   else; cmap = varargin{1}; end
            if nargin<4 || isempty(varargin{2});  figNum = 1;        else; figNum = varargin{2}; end
            if nargin<5 || isempty(varargin{3});  evalstr = '';      else; evalstr = varargin{3}; end

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
    end
end
