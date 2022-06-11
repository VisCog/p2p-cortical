% p2p_psy
%
%  holds all support functions for psychophysical temporal data
%
% functions can be called from outside with 'p2p_psy.<function name>'


classdef p2p_psy
    methods(Static)
        %% pulse width data
        function electrode = Dobelle74_getDataPW()
            pw = [ 1000 500 250 125 62] * 10.^-6; % pulse width (s)
            ind = find(pw ==1000* 10.^-6);
            amp= [4  5 NaN 6 NaN;  2  3 5 8 NaN; 2 3 5 8 NaN; 2 3 5 7 12; 3 3 5 6 NaN]';%milliamps
            for e = 1:size(amp, 2) % cols is the individual electrodes
                for i = 1:length(pw)
                    electrode(e).trl(i).amp = 1000 * amp(i, e);
                    electrode(e).trl(i).amp_n = amp(i, e)./amp(ind, e);
                    electrode(e).trl(i).pw = pw(i);
                    electrode(e).trl(i).freq = 50;
                    electrode(e).trl(i).dur = 1; % don't actually know the duration, not in the paper!
                end
                electrode(e).norm_amp = electrode(e).trl(ind).amp;
                electrode(e).norm_ind = ind;
                electrode(e).color = 'r';
            end
        end
        function electrode = Henderson1979_getDataPW()
            pw = [250 500 750 1000 1500] * 10.^-6; % pulse width (s)
            ind = find(pw ==1000* 10.^-6);
            amp = [4.09 NaN NaN NaN; 2.9 2.8 2.7 2.65;  2.43 2.35 2.29 NaN; 2.1 1.89 NaN NaN; 1.62 1.57 NaN NaN];
            for e = 1:size(amp, 2) % cols is the individual electrodes
                for i = 1:length(pw)
                    electrode(e).trl(i).amp = amp(i, e)*1000;
                    electrode(e).trl(i).amp_n = amp(i, e)./amp(ind, e);
                    electrode(e).trl(i).pw = pw(i);
                    electrode(e).trl(i).freq = 50;
                    electrode(e).trl(i).dur = 1; % don't actually know the duration, not in the paper!
                end
                electrode(e).norm_amp = electrode(1).trl(ind).amp;
                electrode(e).norm_ind = ind;
                electrode(e).color = 'm';
            end
        end

        function electrode = Dobelle79_getDataPW()
            % these data also reported (the 0.5 s pulse train data) in Girvin 1979)
            pw = [125 250 500 1000 2000] * 10.^-6; % pulse width (s)
            amp= [ 2.28  1.73 1.49 1.13 0.77];
            ind = find(pw ==1000* 10.^-6);
            for i = 1:length(pw)
                electrode.trl(i).amp = 1000 * amp(i);
                electrode.trl(i).amp_n = amp(i)./amp(ind);
                electrode.trl(i).pw = pw(i);
                electrode.trl(i).freq = 50;
                electrode.trl(i).dur = .5;
            end
            electrode.norm_amp = electrode.trl(ind).amp;
            electrode.norm_ind = ind;
            electrode.color = 'c';
        end
        function electrode = Brindley1968_getDataPW()
            amp= [8 9 9 10 13 16 19 25 28 36 56]; % threshold volts
            pw = [ 1000 600 400 300 200 100 60 40 30 20 10] * 10.^-6; % pulse width (s)
            ind = find(pw ==1000 * 10.^-6);
            W2A = 1/3;  % estimated resistance of 3000, convert to microamp
            for i = 1:size(amp,2)
                electrode.trl(i).amp = 1000 *  W2A * amp(i);
                electrode.trl(i).amp_n = amp(i)./amp(ind);
                electrode.trl(i).pw = pw(i);
                electrode.trl(i).freq =  30;
                electrode.trl(i).dur = .5; % don't actually know the duration, not in the paper!
            end
            electrode.norm_amp = electrode.trl(ind).amp; % this is the current for a 1000usec pulse, which is our 'standard'
            electrode.norm_ind = ind; % this is the index for the 1000us pulse
            electrode.color = 'b';
        end

        function electrode = Girvin79_getDataPW()
            % these data also reported (the 0.5 s pulse train data) in Girvin 1979)
            pw = [250 500 1000 2000] * 10.^-6; % pulse width (s)
            amp = [3.85 3 2.24 1.84];
            ind = find(pw ==1000* 10.^-6);
            for i = 1:length(pw)
                electrode.trl(i).amp = 1000 * amp(i);
                electrode.trl(i).amp_n = amp(i)./amp(ind);
                electrode.trl(i).pw = pw(i);
                electrode.trl(i).freq = 1;
                electrode.trl(i).dur = .5;
            end
            electrode.norm_amp = electrode.trl(ind).amp;
            electrode.norm_ind = ind;
            electrode.color = 'g';
        end


        function electrode = Fernandez2021_getDataPW()
            pw = [ 100 170 400 800]*10.^-6; %(s)
            charge = [8.065 6.645 16.128 20.4]/1000;
            amp = charge./pw;
            ind = length(pw);
            for i = 1:length(pw)
                electrode.trl(i).amp = amp(i);
                electrode.trl(i).amp_n = amp(i)./amp(ind); % never did a 100ms pulse width so can't normalize
                electrode.trl(i).pw = pw(i);
                electrode.trl(i).freq = 300;
                electrode.trl(i).dur =  0.1666;
                electrode.trl(i).ip = 60*10.^-6;
            end
            electrode.norm_amp = electrode.trl(ind).amp;
            electrode.norm_ind = ind;
            electrode.color = 'k';
            disp('WARNING - Fernandez 2021');
            disp('1000us not available for normalizing, using 800us instead');
        end
        function electrode = Tehovnik2007_getDataPW()
            pw = [0.04 0.05 0.07 0.09 0.13 0.19 0.26 0.35 0.45 0.59 0.68 0.76 0.79 0.79 ] * 10.^-3; % pulse width (s)
            amp = [5.72 4.95 4.06 3.27 2.58 1.97 1.65 1.26 1.03 0.91 0.88 0.84 0.84 0.84];
            ind = length(amp);
            for i = 1:length(pw)
                electrode.trl(i).amp = amp(i) * 1000;
                electrode.trl(i).amp_n = amp(i)./amp(ind); % never did a 100ms pulse width so can't normalize
                electrode.trl(i).pw = pw(i);
            end
            electrode.norm_amp = electrode.trl(ind).amp;
            electrode.norm_ind = ind;
            electrode.color = 'y';
            disp('WARNING - Tehovnik 2007');
            disp('1000us not available for normalizing, using 790us instead');
        end

        %% frequency data
        function electrode = Brindley1968_getDataFreq()
            amp = [29 27 21 21 25 35 37 39 35 29]; % threshold current
            freq = [25 50 100 160 250 400 630 1000 1600 4000]; % freq
            ind = find(freq == 250);
            W2A = 1/3;
            for i = 1:size(amp,2)
                electrode.trl(i).amp = 1000 * amp(i) * W2A; % estimated resistance of 3000
                electrode.trl(i).amp_n = amp(i)./amp(ind);
                electrode.trl(i).pw = 30 * 10.^-6;
                electrode.trl(i).freq = freq(i);
                electrode.trl(i).dur = .5; % probably this, see Girvin
                electrode.trl(i).ip = 60*10.^-6;
            end
            electrode.norm_amp = electrode.trl(ind).amp;
            electrode.norm_ind = ind;
            electrode.color = 'b';
            disp('WARNING - Brindldey 1968');
            disp('200 Hz  not available for normalizing, using 250 Hz instead');
        end

        function electrode = Dobelle74_getDataFreq()
            freq = [200 50 25 12];
            ind = find(freq == 200);
            % cols is the individual electrodes, 10 electrode 4 frequencies
            amp= [5 3 2 2 2 2 2 2 2 2;  5 4 4 4 4 4 2 2 4 2;NaN 4 NaN 4 NaN NaN NaN 3 NaN NaN;
                NaN 8 NaN NaN NaN NaN NaN NaN NaN NaN];
            for e = 1:size(amp, 2)
                for i = 1:length(freq)
                    electrode(e).trl(i).freq = freq(i);
                    electrode(e).trl(i).amp = amp(i,e) * 1000; % estimated resistance of 3000
                    electrode(e).trl(i).amp_n = amp(i,e)./amp(ind,e);
                    electrode(e).trl(i).pw = 500 * 10.^-6;
                    electrode(e).trl(i).dur = 1; % don't actually know the duration, not in the paper!
                end
                electrode(e).norm_amp =  electrode(e).trl(ind).amp;
                electrode(e).norm_ind = ind;
                electrode(e).color = 'r';
            end
        end

        function electrode = Girvin79_getDataFreq()
            % these data also reported (the 0.5 s pulse train data) in Girvin 1979)
            freq = [12.5 25 50 100 200 400 800 1600];
            amp = [2.72 2.20 1.90 1.43 1.04 0.89 0.86 0.81 ];
            ind = find(freq == 200);
            for i = 1:length(freq)
                electrode.trl(i).amp = amp(i) * 1000;
                electrode.trl(i).amp_n = amp(i)./amp(ind);
                electrode.trl(i).pw = 250 * 10.^-6;
                electrode.trl(i).freq = freq(i);
                electrode.trl(i).dur = .5;
            end
            electrode.norm_amp = electrode.trl(ind).amp;
            electrode.norm_ind = ind;
            electrode.color = 'g';
        end
        function electrode = Fernandez2021_getDataFreq()
            % these data also reported (the 0.5 s pulse train data) in Girvin 1979)
            freq = [300 200 100];
            charge =  [6.57 9.39 14.18]./1000;
            amp = charge/(170*10.^-6);
            ind = find(freq == 200);
            for i = 1:length(freq)
                electrode.trl(i).amp = amp(i);
                electrode.trl(i).amp_n = amp(i)./amp(ind);
                electrode.trl(i).pw = 170*10.^-6;
                electrode.trl(i).freq = freq(i);
                electrode.trl(i).dur = .1666;
            end
            electrode.norm_amp = electrode.trl(ind).amp;
            electrode.norm_ind = ind;
            electrode.color = 'k';
        end

        %% full psychometric functions
        function [freq, dur, pw]= Fernandez_Data()
            % figure 2A
            pw{1}.val  = 800; pw{1}.col = 'r';
            pw{2}.val  = 400; pw{2}.col = 'g';
            pw{3}.val  = 170; pw{3}.col = 'm';
            pw{4}.val  = 100; pw{4}.col = 'b';
            pw{1}.name = 'pw'; pw{1}.lim = [0 40];
            pw{1}.data = [0   0.1;   7.790   0.1; 15.980   0.220; 23.930   0.810; 32.05   0.8; 40   0.910];
            pw{2}.data = [ 0.08  0; 4.1   0.09;  7.870  0;  11.97   0.290;  16.15   0.5; ...
                20.160   0.690;  23.930   1.0;  28.0300   1.000;  32.05   0.9;  35.98   1.0; 40.0   1.0];
            pw{3}.data = [0   0;   1.560  0;  3.440   0.090;   5.0   0.3;   6.8   0.5; ...
                8.440   0.91;  10.160   0.89;  11.97   0.8;  13.690   1.0;  15.490   1.01;  16.89   1.0];
            pw{4}.data = [0.9  0;  1.89   0.1;   2.79 0;   3.93 0;  4.920   0.2
                6.070   0.3;  6.890   0.3;  8   0.6;  9  0.710;   10  0.8];

            % figure 2B
            freq{1}.val  = 100; freq{1}.col = 'b';
            freq{2}.val  = 200; freq{2}.col = 'g';
            freq{3}.val  = 300; freq{3}.col = 'r';
            freq{1}.name = 'freq'; freq{1}.lim = [0 16];
            % charge per phase, propn seen
            freq{1}.data = [ 3.41  0;   5.12    0;  6.83   0; 8.38  0.4; ...
                10.22 0.29;  11.84  0.39;  13.5  0.4; 16.87   0.49];
            freq{2}.data = [  0   0.20;  1.78   0.19; 3.37 0;  5.16 0; ...
                6.78  0.09;  8.46  0.6;  11.96   0.59;  13.4   0.65;  15.23   0.74];
            freq{3}.data = [0  0;  1.71  0;   3.36   0.09; 5.02   0.3000;   6.71  0.50; ...
                8.4   0.91;  10.19  0.89; 11.86 0.80;13.51  1.00;15.18  0.99; 16.88   0.99];

            dur{1}.val  = 50; dur{1}.col = 'b';
            dur{2}.val  = 83; dur{2}.col = 'g';
            dur{3}.val  = 166; dur{3}.col = 'r';
            dur{1}.name = 'dur'; dur{1}.lim = [0 16];
            dur{1}.data = [0   0;   1.480  0;  3.13  0;  4.92 0;  6.49  0.12; ...
                8.33  0.11;   9.95   0.22; 11.70   0.22;  13.48   0.51;  15.16   0.78; 16.84  0.61];
            dur{2}.data=[  0   0.21;  3.25  0.21;4.83   0.01;  6.62   0.40;  8.34   0.56; ...
                10.15   0.70; 11.80   0.70;  15.21   0.80; 16.91   0.80];
            dur{3}.data=[ 0   0.210;    1.48   0;  3.25   0.2;  4.82   0.32;  6.54   0.51;  8.38   0.9; ...
                10.18   0.88;  13.35   1.0;  15.19   1.0; 16.94   1.0];

            brightness.data =  [10   0.00140 0.0003  0.0024*1.0e+02;
                20   0.001700 0.0004 0.00290*1.0e+02; ...
                30  0.004800  0.0016 0.0084*1.0e+02; ...
                40   0.008600  0.0042 0.0136*1.0e+02; ...
                50   0.012   0.0082 0.0160*1.0e+02; ...
                60   0.020200  0.017 0.0231*1.0e+02; ...
                70  0.020500   0.01790 0.0233*1.0e+02; ...
                80   0.020300   0.0166 0.02410*1.0e+02; ...
                90   0.025500  0.0225 0.02830*1.0e+02; ...
                100   0.023   0.0197 0.0261];
        end
    end
end
