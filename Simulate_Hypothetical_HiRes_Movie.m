% Simulate_Hypothetical_HiRes_Movie.m
% Simulates what a movie would look like using an
% idealized high resolution array
%
% written IF & GMB
% 25/02/2023 moved into clean folder (IF)

clear
n_electrodes = 4380;
filename = 'CatNotCooperating';
tp = p2p_c.define_temporalparameters();

%% first quickly check we can load all the receptive fields mat files
if 0
    for e = 1:4255
        filename = ['datasets'  filesep 'Optimal_Spacing' filesep 'Optimal_Spacing_RFmaps_bosking_', num2str(e), '.mat'];
        try
            tmp = load(filename);
        catch
            disp(['file ', num2str(e), ' missing']);
        end
    end
end

%% load input movie
filename_in =['movies' filesep filename '.avi'];
vid_in = VideoReader(filename_in);
vid_dur  = vid_in.NumFrames;

%% open output video file
filename_out = ['movies' filesep filename '_Optimal_Spacing_hires.avi'];
vid_out = VideoWriter([filename_out, '.avi']);
vid_out.FrameRate = vid_in.FrameRate;
open(vid_out);

%% load saved rfmaps
filename = ['datasets'  filesep 'Optimal_Spacing' filesep 'Optimal_Spacing_RFmaps_bosking_1.mat'];
tmp = load(filename);
rfmap = tmp.rfmap;

for k = 1:vid_dur% for each frame

    disp(['simulating frame ', num2str(k), ' out of ', num2str(vid_dur)]);
    frame = readFrame(vid_in);
    size(frame, 1)~=size(frame, 2)
    frame = frame(:, 421:1500, :); 
    if k==20
        frame = mean(frame, 3);
        in_movie=imresize(frame, [size(rfmap, 1), size(rfmap, 2)], "bilinear");
        in_movie  = in_movie./255; % scale the movie between 0 -1
        img = zeros(size(rfmap)); B0 = zeros(size(rfmap));
        ct = 0;
        for e = 1:1:n_electrodes
            if mod(e, 100)==1;  disp(['simulating ', num2str(e), ' out of ', num2str(n_electrodes)]);   end
            filename = ['datasets'  filesep 'Optimal_Spacing' filesep 'Optimal_Spacing_RFmaps_bosking_', num2str(e), '.mat'];
            tmp = load(filename);
            rfmap = tmp.rfmap./sum(tmp.rfmap(:)); %each rf normalized to have an area of 1
            if ~isnan(rfmap)
                amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
                %     tp.model = 'compression';
                img = img + (rfmap.*amp);
                B0 =  B0 + (rfmap*1); % deals with non-homogeneity in the distribution of rfs
                ct = ct+1;
            end

        end
        disp(ct)
        % weirdly, because we're not really simulating the temporal part of the
        % model that relates electrical stimulation to brightness we don't include the nonlinearity either
        % we assume we want the brightness to be related to the image intensity
        figure(2)
        B0_norm = B0;
        img_norm =  img;
        img_norm(img_norm<=0 | B0_norm<=0) = 0;
        B0_norm(B0_norm<=0 | img_norm<=0)=0;
        norm = img_norm./B0_norm;
        norm(norm>1) = 1;
        subplot(1,2,1)
        imagesc(in_movie); colormap(gray(256));
        axis equal;    axis tight;    axis off
        subplot(1,2,2)
        image(norm*255); colormap(gray(256));
        axis equal; axis off;  axis tight
        frame = getframe(gca);
        writeVideo(vid_out,frame);
    end
end
close(vid_out);
