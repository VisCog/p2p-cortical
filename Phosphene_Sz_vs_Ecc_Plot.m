% Simulate_Sz_vs_Ecc_ESize
% plots data examining phosphene size as a function of
% electrode size and eccentricity
%
% written IF 05/03/2023
clear 
figure(1);clf

T = readtable('datasets/PhospheneSz_vs_Ecc_keliris_newsize.xlsx');
subplot(1,2,1)
szList = unique(T.esize);
eccList = unique(T.eccentricity);
clist = redgreen(length(szList));
for s = 1:length(szList)
    for ecc = 1:length(eccList)
        eid =T.esize==szList(s) ;
        Tid = T(eid, :);

    end
            h(s)=  plot(Tid.eccentricity, Tid.sz_ellipse*2, '.-',  ....
           'Color', clist(s,:), 'MarkerSize', 18); hold on
    if s==1
        linfit_keliris = polyfit(Tid.eccentricity, Tid.sz_sigma, 1); % use the smaller sigma size for packing
    end
end
title('Keliris'); xlabel('eccentricity');
legend(h, num2str(round(szList/1000, 2)))
axis([ 0 37 0 22])


% so the simulation gets a bit blooey near the fovea (takes infinitely long
% to do simulations at that resolution


subplot(1,2,2)
T = readtable('datasets/PhospheneSz_vs_Ecc_bosking_newsize_si.xlsx');
szList = unique(T.esize);
eccList = unique(T.eccentricity);
clist = redgreen(length(szList));
for s = 1:length(szList)
    for ecc = 1:length(eccList)
        eid =T.esize==szList(s) ;
        Tid = T(eid, :);

    end
            h(s)=  plot(Tid.eccentricity, Tid.sz_ellipse*2, '.-',  ....
           'Color', clist(s,:), 'MarkerSize', 18); hold on
        if s==1
        linfit_bosking = polyfit(Tid.eccentricity, Tid.sz_sigma, 1); % use the smaller sigma size for packing
    end
end
plot(Tid.eccentricity, (0.1787 +Tid.eccentricity*0.262) , '--','Color', 'b', 'LineWidth', 1)
title('Bosking'); xlabel('eccentricity');
axis([ 0 37 0 22])

