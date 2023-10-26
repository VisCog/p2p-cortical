
function array = ArraySim_main(array)
% Calls the functions needed to simulate various arrays and examine what
% the perceptual experience would be.\
% send in a, must contain:
% array.spaceFac = 4 (or 3, 2, 1)
% array.arrayList =  {'optimal', ' regular_cortex','regular_visualfield'};%


if strcmp(computer, 'PCWIN64')
    params.homedir = ('C:\Users\Ione Fine\Documents\code\p2p-cortical_new');
    params.datadir = [params.homedir,   filesep, 'datasets', filesep, 'ArrayRFmaps'];
    par_size=1;
elseif strcmp(computer,  'GLNXA64') && array.isbeeble
    params.homedir =  '/home/ionefine/Documents/code/p2p-cortical_new';
    params.datadir = [params.homedir,   filesep, 'datasets', filesep, 'ArrayRFmaps'];
    par_size = 6;
elseif  strcmp(computer,  'GLNXA64') && ~array.isbeeble
    params.homedir = ['/mnt', filesep, 'viscog', filesep, 'FineLab', filesep, 'IoneFine', filesep, 'code', filesep, 'p2p-cortical_new'];
    params.datadir = ['/mnt', filesep, 'viscog', filesep, 'FineLab', filesep, 'IoneFine', filesep, 'code', filesep, 'p2p-cortical_new', ...
        filesep, 'datasets', filesep, 'ArrayRFmaps'];
    par_size = 12;
end
disp(['using a parallel pool of ', num2str(par_size)]);
if array.spaceFac == 4
    array.nelect = 79;
    c.pixpermm = 16; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =16; %16; % 18;
elseif  array.spaceFac ==3;
    array.nelect = 163;
    c.pixpermm = 16; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =16; %16; % 18;
elseif array.spaceFac ==2;
    array.nelect = 399;
    c.pixpermm = 14; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =14; %16; % 18;
elseif array.spaceFac == 1
    array.nelect = 1880;
    c.pixpermm = 10; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =12; %16; % 18;
end

array.eccLim = [0  32];
c.cortexHeight = [-40,40]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 80];
if array.spaceFac <3

else
    c.pixpermm = 16; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =16; %16; % 18;
end
if strcmp(array.esize, 'large')
    disp('simulating large electrodes')
    c.I_k =6.75;
    c.radius = 0.25;
else
    c.I_k = 1000;
 c.radius = 0.001;
end
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-array.eccLim(2)*1.5,array.eccLim(2)*1.5];
v.visfieldWidth=[-array.eccLim(2)*1.5,array.eccLim(2)*1.5];


params.plot = 1;
params.flag_flip = [1, 1 ]; % up down then left right
params.scale = [1 127]; % converts to uint8
params.overwrite = 1;
arrayList = array.arrayList;
if par_size>1
    delete(gcp)
    parpool(par_size);
end

for aa = 1:length(array.arrayList)
    array.arrayStyle = arrayList{aa}; % hack so can't get overwritten
    array.filename = [params.datadir, filesep, array.arrayStyle, filesep, 'Array_Sim_', array.arrayStyle, '_', num2str(array.spaceFac), '_', array.esize];
   
    %% calculate locations
    if 1
        array = p2p_c.Array_Sim_Location(c, array);
        save([array.filename,'.mat'], 'array','-mat');
    else % load a previously generated file
        array = load([array.filename,'.mat'], '-.mat');
    end

    %% plot locations
    if 1
        p.plotNums = [1, 3, 4, 2];
        p2p_c.Array_Sim_Plot(array, c, p);
         figure(3)
        set(gca, 'XLim', [-48 48]);
         set(gca, 'YLim', [-48 48]);
                 figure(4)
         set(gca, 'XLim', [-5 80]);
         set(gca, 'YLim', [-30 30]);
    end

    %% generate RFs
    if 0
        array.T = array.Tfull;
        % if we are flipping the array then we may not compute for all the
        % electrodes
        if params.flag_flip(1)==1;  array.T  = array.T (array.T.vy>0, :); end %up down
        if params.flag_flip(2)==1;  array.T  = array.T (array.T.vx<0, :); % left right
        end

        tmp = par_size*ceil(height(array.T)/par_size);
        range = reshape(1:tmp, tmp/par_size, par_size);
        if par_size==1
            for r = 1:par_size
                rr = range(:, r);
                rfmaps = p2p_c.Array_Sim_GenerateRFmaps(array, v, c,params, rr);
                subset(r).rfmaps = rfmaps;
            end
        else
            parfor r = 1:par_size
                rr = range(:, r);
                rfmaps = p2p_c.Array_Sim_GenerateRFmaps(array, v, c,params, rr);
                subset(r).rfmaps = rfmaps;
            end
        end

        array.rfmaps = subset(1).rfmaps;
        for r = 2:length(subset)
            array.rfmaps = cat(3, array.rfmaps, subset(r).rfmaps);
        end
        disp(['Finished RF creation ', array.arrayStyle]);
        save([array.filename, '_rfs.mat'], 'array', 'params', '-v7.3')
    end
    % put the maps back into the matrix
end
end

