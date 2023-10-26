function ArraySimulation_main_movie(a)
% Calls the functions needed to simulate various arrays and examine what
% the perceptual experience would be.
tic
if strcmp(computer, 'PCWIN64')
    params.homedir = ('C:\Users\Ione Fine\Documents\code\p2p-cortical_new');
    params.datadir = [params.homedir,   filesep, 'datasets', filesep, 'ArrayRFmaps'];
    %
    %     params.datadir = ['Z:', filesep, 'viscog', filesep, 'FineLab', filesep, 'IoneFine', filesep, 'code', filesep, 'p2p-cortical_new', ...
    %     filesep, 'datasets', filesep, 'ArrayRFmaps',filesep];
    par_size = 1;
else
    params.homedir = ['/mnt', filesep, 'viscog', filesep, 'FineLab', filesep, 'IoneFine', filesep, 'code', filesep, 'p2p-cortical_new'];
    params.datadir = ['/mnt', filesep, 'viscog', filesep, 'FineLab', filesep, 'IoneFine', filesep, 'code', filesep, 'p2p-cortical_new', ...
        filesep, 'datasets', filesep, 'ArrayRFmaps'];
    par_size = 10;
    parpool(par_size);
end

a.filename = [params.datadir, filesep, a.arrayStyle, filesep, 'Array_Sim_loc_', a.arrayStyle, '_', num2str(a.spaceFac)];

load([a.filename, '_rfs', '.mat'])
% for r = 1:length(subset)
%     tmp = subset(r).a.rfmaps;
%     a.rfmaps(:, :, range(:, r)) = tmp;
% end
%% create movie

m.keepframes = [1, Inf]; params.flag_flip = [1, 1];
m.crop = 70;
p2p_c.Array_Spacing_Movie(a,params, m)
end

