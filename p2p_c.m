% p2p_c
%
%  holds all support functions for p2p cortex
% project.
%
% functions can be called from outside with 'p2p_c.<function name>'


classdef p2p_c
    methods(Static)
        %% definitions - temporal properties
        function tp = define_temporalparameters(varargin)
            % sets up the parameters for the temporal model

            if nargin>0
                tp = varargin{1};
            else
                tp = [];
            end
            if ~isfield(tp, 'dt'); tp.dt = .001 * 10^-3; end % 001 time sampling in ms, should be no larger than 1/10 of tau1
            if ~isfield(tp, 'tau1'); tp.tau1 = 0.178 * 10^-3; end
            % .0.12 * 10^-3 fit based on Brindley, 1968
            % Tehovnk 2004 estimates 0.13-0.24 ms
            % Fernandez 2021 looks like 0.178
            if ~isfield(tp, 'tau2'); tp.tau2 =  26.250* 10^-3; end % 24-33 from retina
            if ~isfield(tp, 'ncascades'); tp.ncascades = 3; end
            % number of cascades in the slow filter of the temporal convolution
            if ~isfield(tp, 'scaleR1'); tp.scaleR1 = 1; end % scaling of current on the retina

            %  the units are current on the retina.
            % Based on Fernandez, for an 800 microsec pulse train at 300Hz the threshold current is 25 mAmps
            % For a 170 microsec pulse at 300 Hz, threshold is 40, which
            % gives 25 at the retina, and saturation is at ~90 (144 at electrode)

            % nonlinearity parameters
            if ~isfield(tp, 'model');  tp.model = 'weibull'; end

            if strcmp(tp.model, 'weibull')
                disp('using weibull semisaturation constant')
                tp.asymptote = 2.5;
                tp.a = 50;
                tp.b = 4;
                tp.scaleR4 = 1; % amount to scale the response by
            elseif strcmp(tp.model, 'sigmoid')
                disp('using sigmoid semisaturation constant')
                tp.asymptote = 2.5;
                tp.sigma = 30; % electrical semisaturation constant
                tp.scaleR4 = 1; % amount to scale the response by
            elseif strcmp(tp.model, 'exp')
                disp('using exponential semisaturation constant')
                tp.asymptote = 2.5;
                tp.k = 0.02;
                tp.scaleR4 = 1;
            elseif strcmp(tp.model, 'linear')
                disp('no semisaturation constant')
                tp.scaleR4 = 1;
            elseif strcmp(tp.model, 'leakyintegrator')
                disp('leaky integrator, no semisaturation constant or slow integration')
                tp.scaleR4 = 1;
            end
        end
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
            if ~isfield(trl, 'lag');    trl.lag = trl.pw;       end % delay before the pulse train begins
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
                lag = round(trl.lag/tp.dt); % delay before beginning the pulse train
                on =  mod(trl.t,1/trl.freq) < trl.pw;
                delay =  trl.pw+trl.ip; % time difference between on and off
                off = mod(trl.t-delay,1/trl.freq) < trl.pw;
                tmp  = trl.amp.*(on-off);
                trl.pt= zeros(1, lag+length(tmp));
                trl.pt(lag+1:lag+length(tmp))=tmp;
            end
            if trl.dur<trl.simdur % usually we simulate a little longer than the trial, to allow for the response
                trl.pt((end+1):round((trl.simdur/tp.dt))) = 0;
                trl.t = 0:tp.dt:trl.simdur-tp.dt;
            end
        end
        function trl = finite_element(tp, trl)
            % Implements accumulation of current over time using a simple finite element method
            % written GMB 11/10/2017
            % adapted for cortex IF 3/2/2018
            % Finite difference method:

            %  tmp.ca = 0;
            %  tmp.cl = 0;
            tmp.R1 = 0;
            tmp.R4 =  zeros(1,tp.ncascades+1);

            for i=1:length(trl.pt)
                tmp.R1 = tmp.R1 + tp.dt * ((tp.scaleR1.*trl.pt(i))-tmp.R1)/tp.tau1;
                % except for very short pulses max R2 is roughly equal to
                % input current
                trl.R1(i) = tmp.R1;
                trl.R2(i) = max(tmp.R1, 0);
                if strcmp(tp.model, 'sigmoid')
                    tmp.R3 = tp.asymptote .* tmp.R2.^2./(tmp.R2.^2 + tp.sigma.^2);
                elseif strcmp(tp.model, 'normcdf')
                    tmp.R3 = tp.asymptote.*normcdf(tmp.R2, tp.mean, tp.sigma);
                    % tmp.R3 = tmp.R3 - tp.asymptote.*normcdf(0, tp.mean, tp.sigma);
                    tmp.R3(tmp.R3<0) = 0;
                elseif strcmp(tp.model,'weibull')
                    tmp.R3 = p2p_c.Weibull(tp,tmp.R2);
                elseif strcmp(tp.model,'linear') 
                    tmp.R3 = tmp.R2;
                end
                if ~strcmp(tp.model, 'simpleleakyintegrator')
                    trl.R3(i) = tmp.R3;
                    tmp.R4(:, 1) = tmp.R3;
                    for j=1:tp.ncascades
                        tmp.R4(:, j+1) = tmp.R4(:, j+1) + tp.dt*(tmp.R4(:, j) - tmp.R4(:, j+1))/tp.tau2;
                    end
                    trl.resp(i)= tp.scaleR4 * tmp.R4(:, end);
                else % if just a simple leaky integrator
                    trl.resp(i) = tp.scaleR4 * trl.R2(i);
                end
            end
        end

        %% definitions, psychophysics
        function amp = find_threshold(trl, tp, varargin)
            % takes in trial, tp, and optional fitParams
            % 
            % finds and returns the trl.amp for which the max output of the model for that trial, trial.resp, is equal
            % to fitParams.thr
            
            if nargin<3; fitParams.tol = 0.001;   fitParams.lo = 0;  fitParams.hi = trl.amp * 10;  fitParams.thr = 1;   fitParams.nreps = 10;
            else;       fitParams = varargin{1};          end

            tol = fitParams.tol; lo = fitParams.lo; hi = fitParams.hi;

            for i = 1:fitParams.nreps
                trl.amp  = (hi+lo)/2;
                trl = p2p_c.define_trial(tp,trl);
                trl = p2p_c.finite_element(tp, trl);

                if max(trl.resp(:)) > fitParams.thr
                    hi = trl.amp;
                else
                    lo = trl.amp;
                end

                if abs(max(trl.resp(:))-fitParams.thr) < tol % close enough
                    disp(['trl.amp  = ', num2str((hi+lo)/2), ' fitParams.thr = ', num2str(fitParams.thr), ' trl.resp = ', num2str(max(trl.resp))]);
                    amp = (hi+lo)/2;
                    disp('close enough ...'); break;
                end
                disp(['trl.amp  = ', num2str((hi+lo)/2), ' fitParams.thr = ', num2str(fitParams.thr), ' trl.resp = ', num2str(max(trl.resp))]);

                % if the fitting doesn't seem to work, check that the original hi and
                % lo ranges aren't too narrow
            end
            amp = (hi+lo)/2;
        end

    end
end