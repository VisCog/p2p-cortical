clear
close all

%% getting movies from cat videos
fname{1}.name  ='C:\Users\Ione Fine\Documents\code\p2p-cortical_new\movies\CatNotCooperating.avi';
fname{2}.name  ='C:\Users\Ione Fine\Documents\code\p2p-cortical_new\datasets\ArrayRFmaps\regular_visualfield\Array_Sim_regular_visualfield_2_tiny_CatNotCooperating.avi';
fname{3}.name  ='C:\Users\Ione Fine\Documents\code\p2p-cortical_new\datasets\ArrayRFmaps\regular_cortex\Array_Sim_regular_cortex_2_tiny_CatNotCooperating.avi';
fname{4}.name  ='C:\Users\Ione Fine\Documents\code\p2p-cortical_new\datasets\ArrayRFmaps\optimal\Array_Sim_optimal_2_tiny_CatNotCooperating.avi';

%fname{1}.name  ='C:\Users\Ione Fine\Documents\code\p2p-cortical_new\datasets\ArrayRFmaps\optimal\Array_Sim_optimal_1_tiny_CatNotCooperating.avi';
%fname{2}.name  ='C:\Users\Ione Fine\Documents\code\p2p-cortical_new\datasets\ArrayRFmaps\optimal\Array_Sim_optimal_2_tiny_CatNotCooperating.avi';

% for movies of simulations
for i = 1:length(fname)
    v = VideoReader(fname{i}.name);
    v.CurrentTime = 0.6;
    vidFrame = readFrame(v);
    subplot(1, 4, i)
    if ndims(vidFrame)>2
        vidFrame = mean(vidFrame, 3);
        if size(vidFrame, 1)<size(vidFrame, 2) % if the image isn't square
            off = 1+ ( (size(vidFrame, 2)-size(vidFrame, 1))/2);
            vidFrame = vidFrame(:, off-1:size(vidFrame, 2)-off, :);
        end
    end
    colormap("gray");
    imagesc(vidFrame);
    axis off
    axis equal
    pause
end

return
%% getting frames from letter movies

fname{1}.name  = 'C:\Users\Ione Fine\Documents\code\p2p-cortical_new\movies\Beauchamp_Cell2020_Figure 4-gridCSimultaneous.avi';
fname{2}.name  = 'C:\Users\Ione Fine\Documents\code\p2p-cortical_new\movies\Beauchamp_Cell2020_Figure 4-gridNSimultaneous.avi';
fname{3}.name  = 'C:\Users\Ione Fine\Documents\code\p2p-cortical_new\movies\Beauchamp_Cell2020_Figure 4-gridSSimultaneous.avi';
fname{4}.name  = 'C:\Users\Ione Fine\Documents\code\p2p-cortical_new\movies\Beauchamp_Cell2020_Figure 4-gridUSimultaneous.avi';
for i = 1:length(fname)
    v = VideoReader(fname{i}.name);
   % while(hasFrame(v))
v.CurrentTime = 0.85;
    vidFrame = readFrame(v);
        imagesc(vidFrame(150:356, 1:end-150, :));      
    axis off
    axis equal; drawnow
    pause
%    end 
end