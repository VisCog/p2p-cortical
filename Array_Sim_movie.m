function array = ArraySim_movie(a)

if  ~isfield(a, 'isbeeble')
    a.isbeeble = 0;
end

if strcmp(computer, 'PCWIN64')
    params.homedir = ('C:\Users\Ione Fine\Documents\code\p2p-cortical_new');
    params.datadir = [params.homedir,   filesep, 'datasets', filesep, 'ArrayRFmaps'];
elseif strcmp(computer,  'GLNXA64') && a.isbeeble
    params.homedir =  '/home/ionefine/Documents/code/p2p-cortical_new';
    params.datadir = [params.homedir,   filesep, 'datasets', filesep, 'ArrayRFmaps'];
elseif  strcmp(computer,  'GLNXA64') && ~array.isbeeble
    params.homedir = ['/mnt', filesep, 'viscog', filesep, 'FineLab', filesep, 'IoneFine', filesep, 'code', filesep, 'p2p-cortical_new'];
    params.datadir = ['/mnt', filesep, 'viscog', filesep, 'FineLab', filesep, 'IoneFine', filesep, 'code', filesep, 'p2p-cortical_new', ...
        filesep, 'datasets', filesep, 'ArrayRFmaps'];
end

params.flag_flip = [1, 1];
m.keepframes = [1, Inf];
m.crop = 70;
arrayList = a.arrayList;

for aa = 1:length(arrayList)
    arrayStyle = arrayList{aa};
    a_filename = [params.datadir, filesep, arrayStyle, filesep, 'Array_Sim_', arrayStyle, '_', num2str(a.spaceFac), '_', a.esize];
    tmp = load([a_filename, '_rfs.mat']);
    array.rfmaps = tmp.array.rfmaps;
    array.spaceFac = tmp.array.spaceFac;
    if array.spaceFac == 4
        array.nelect = 79;
        c.pixpermm = 16; %16; % 16; % creates a rfmap of 936
        params.pixperdeg =16; %16; % 18;
    elseif  array.spaceFac ==3
        array.nelect = 163;
        c.pixpermm = 16; %16; % 16; % creates a rfmap of 936
        params.pixperdeg =16; %16; % 18;
    elseif array.spaceFac ==2
        array.nelect = 399;
        c.pixpermm = 14; %16; % 16; % creates a rfmap of 936
        params.pixperdeg =14; %16; % 18;
    elseif array.spaceFac == 1
        array.nelect = 1880;
        c.pixpermm = 10; %16; % 16; % creates a rfmap of 936
        params.pixperdeg =10; %16; % 18;
    end

    gvals  = ones(1, size(array.rfmaps, 3));
    for j = 1:size(array.rfmaps, 3)
        r = array.rfmaps(:, :, j);
        if max(r)==min(r)
            gvals(j) = 0;
        end
    end
    array.rfmaps =   array.rfmaps(:, :, find(gvals));

    disp(['number of electrodes that produced phosphenes = ', num2str(length(find(gvals)))]);
    params.scale  = tmp.params.scale;

    params.offset = 75; % shift luminance up in movie

    m.filename_in =['movies' filesep 'original_movies', filesep,  'CatNotCooperating.avi'];
    m.filename_out =[a_filename,  '_CatNotCooperating.avi'];
    p2p_c.Array_Sim_Movie(array,params, m);
    disp(['completed ',     arrayStyle, '  ',num2str(a.spaceFac)])
end
end
