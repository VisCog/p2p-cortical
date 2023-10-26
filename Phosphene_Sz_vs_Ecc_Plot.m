% Simulate_Sz_vs_Ecc_ESize
% plots data examining phosphene size as a function of
% electrode size and eccentricity
%
% written IF 05/03/2023
clear
figure(1);clf

T = readtable('PhospheneSz_vs_Ecc_keliris.xlsx');
szList = unique(T.esize);
eccList = unique(T.eccentricity); 
eccList = eccList(find(eccList<34));
clist = viridis(length(szList)+1);
for s = 1:length(szList)
    eid1 = T.esize==szList(s);
    for ecc = 1:length(eccList)
        eid2 = T.eccentricity==eccList(ecc) ;
        Tid = T(find(eid1.*eid2), :);
        if ~isempty(Tid)
            mn_sz(s,ecc) = mean(Tid.sz_ellipse*2);
            std(:, s,ecc) = abs(prctile(Tid.sz_ellipse*2, [5, 95])-mn_sz(s, ecc));
        end
    end
end

for s = 1:length(szList)
    hp=  plot(eccList, mn_sz(s, :), '.-',  ....
        'MarkerFaceColor', clist(s,:), 'MarkerSize', 30, ...
        'LineStyle', '-', 'LineWidth', 2, 'MarkerEdgeColor', clist(s,:), 'Color', clist(s,:)); hold on
    sh =  shadedErrorBar(eccList, mn_sz(s, :),  squeeze(std(:, s, :)));
    sh.patch.FaceColor = clist(s, :);
        sh.patch.FaceColor = clist(s, :);
         sh.mainLine.Color = clist(s, :);
           sh.edge(1).Color ='none';  sh.edge(2).Color = 'none';
    if s==1
        linfit_keliris = polyfit(Tid.eccentricity, Tid.sz_sigma, 1); % use the smaller sigma size for packing
    end
end
title('Keliris'); xlabel('eccentricity');
legend(hp, num2str(round(szList/1000, 2)))
axis([ 0 27 0 22])

