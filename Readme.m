% readme

%    p2p_c.m : the guts of the program, a collection of methods
  

%% spatial stuff
% 
% SimulateMaps.m  
% creates figure showing orientation, ocular dominance maps etc

%  Winawer_PhosphenePictures_Neuron2016.m   
%  creates pictures of phosphenes designed to match
% Winawer Neuron 2016 paper. Uses ugly code to replot Winawer (I'd
% recommend going directly to their code which is publically available if
% you want to plot their data).

% OptimalSpacing
% OptimalSpacing_main - calls the various sub programs
% Optimal Spacing_Location - finds locations for optimal placement of
% electrodes, can vary the spacing
% OptimalSpacing_RFmap - generates the maps for those locations
% OptimalSpacing_Movie - generates a cat movie, based on RF maps

%% temporal maps
% Plot_Temporal_Model.m
% creates figure describing the temporal model       

% Fit_Frequency_Curves.m
% threshold as a function of frequency

% Fit_PulseWidth_Curves.m
% creates chronaxie figure, threshold as a
% function of pulse width

% Winawer_Brightness_Neuron2016.m    
% plots phosphene brightness as a function of pulse parameters

%%  spatiotemporal stuff
% Bosking_SizeAmplitude_JNeuro2017.m 
% plots phosphene size as function of amplitude

% Beauchamp_DynamicLetters_Cell2020_RF.m
% generates the receptive field maps for each of Beauchamps electrodes,
% Figure 4 & 6

% Beauchamp_DynamicLetters_Cell2020_Movie.m
% uses the receptive field maps to generate a movie of letter stimulation

% Winawer_Size_Neuron2016.m   
% plots phosphene size as a function of charge. 




%%  TO DO

%  Simulate_Hypothetical_HiRes_Movie.m       
% Simulate_WFMA.m                           
      
% Bosking_SizeEccentricity_JNeuro2017.m                                      
% Phosphene_Sz_2D.m                                                    
% Phosphene_Sz_vs_Ecc.m                                                
                              
%% datasets
%
% Bosking2017_data
% used grabdata to find size as a function of eccentricity for Bosking et
% al., 2017
                 
              