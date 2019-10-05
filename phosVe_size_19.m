% Bosking_ARVO_19.m

cd('C:\Users\Ione Fine\Documents\code\p2p-cortical\Figures')
rng(1)
c.efthr = 0.01;

v.drawthr = .015;

% Set  electrode locations

eLoc = [.2 .5 .7 1 1.5 2 5 7 10 15 20];
eSize = [.01 .05 .1 .2 .3 .5 1 2 4 8 10];
ct = 1;
for jj=1:length(eLoc)
    for ii=1:length(eSize)
        v.e(ct).ang = 0;
        v.e(ct).ecc = eLoc(jj);
        c.e(ct).radius = eSize(ii);
        ct = ct+1;
    end
end

tp.scFac = 1/10;  %1/10
tp = p2p_c.define_temporalparameters(tp);

sim_sizes = [];
ct = 1;

for jj=1:length(eLoc)
    skip = 0;
    if eLoc(jj)<=1
        c.cortexCenter = [0,0]; c.cortexSize = [20,20]; c.pixpermm = 32;
        v.retinaSize = [6,6];v.pixperdeg = 15;
        
    elseif eLoc(jj)>1 && eLoc(jj)<=6
        c.cortexCenter = [30,0]; c.cortexSize = [60,40]; c.pixpermm = 18;
        v.retinaSize = [20,20];v.pixperdeg = 15;
        
    elseif eLoc(jj)>6
        c.cortexCenter = [30,0]; c.cortexSize = [60,80]; c.pixpermm = 5;
        v.retinaSize = [80,80];v.pixperdeg = 5;
    end
    
    c = p2p_c.define_cortex(c);
    v = p2p_c.define_visualmap(v);
    [c, v] = p2p_c.generate_corticalmap(c, v);
    
    for ii=1:length(eSize)
        if ~(eLoc(jj)<=1 && eSize(ii)>5)
            
        disp(sprintf('Electrode %d of %d',ct,length(v.e)));
        c = p2p_c.define_electrodes(c, v);
        c = p2p_c.generate_ef(c);
        v = p2p_c.generate_rfmap(c, v, ct);
        
        tmp = [];
        tmp.expname = 'generic';
        tmp.amp = 2000; tmp.e = ct;
        tmp = p2p_c.define_trial(tp,tmp);
        
        tmp = p2p_c.generate_phosphene(v, tp, tmp);
        tmp.sim_radius= mean([tmp.ellipse(1).sigma_x tmp.ellipse(1).sigma_y]);
        tmp.sim_diameter = 2 * tmp.sim_radius;
        tmp.sim_brightness = max(tmp.maxphos(:));
        tmp.maxresp = max(tmp.resp);
        p2p_c.plotretgrid(tmp.maxphos(:, :, 1)*1000, v, gray(256), 1,['subplot(2,1,1); title(''phosphene'')';]);
        p2p_c.plotcortgrid(c.e(ct).ef * 256, c, gray(256), 1,['subplot(2,1,2); title(''electric field'')']);
        saveas(gcf, ['corticalFig_eSize', num2str(eSize(ii)*100), '-eLoc_', num2str(eLoc(jj)*100)], 'jpeg');
       
        sim.diameter(ii,jj) =  tmp.sim_diameter;
        sim.loc(ii, jj) =  eLoc(jj);
        sim.esize(ii,jj) =  eSize(ii);
        else
                    sim.diameter(ii,jj) =  NaN;
        sim.loc(ii, jj) =  eLoc(jj);
        sim.esize(ii,jj) =  eSize(ii);
        end
        ct = ct+1;
    end
end


%%
% Plot phosphene size as a function of phosphene eccentricity (for
% stimulation at 1000 microamps)

figure(3)
clf
cmap = hsv(12);
%surf(sim.loc , sim.esize, sim.diameter)
for i = 1:size(sim.esize,1)
    gv = find(~isnan(sim.diameter(i,:)));
    h(i) = plot(sim.loc(i,gv), sim.diameter(i,gv), 'ko-','Color', cmap(i,:),'MarkerFaceColor', cmap(i,:), 'MarkerSize', 8,'LineWidth',1); hold on
end
plot([0 20], [min(sim.diameter(:)), min(sim.diameter(:))], 'k--')
text(12, min(sim.diameter(:))*2, num2str([min(sim.diameter(:))]))
    set(gca, 'XLim', [ 0 21])
    xlabel('eccentricity (deg)')
    ylabel('phosphene diameter (deg)')
    legend(h, eSize)

return

hold on
for jj=1:length(eLoc)
    if eLoc(jj)<=1
        eSize = [.01  .3  .5 1 2];
    elseif eLoc(jj)>1 && eLoc(jj)<=6
        eSize = [.01  .3  .5 1 24 8 10 ];
    elseif eLoc(jj)>6
        eSize = [.01  .3  .5 1 24 8 10 20];
    end
    text(17, sim_sizes(jj, end-1), num2str([eLoc(jj)]));
    
    plot(eSize, sim_sizes(jj,1:length(eSize)), 'ko-','Color', cmap(jj,:),'MarkerFaceColor', cmap(jj,:), 'MarkerSize', 18,'LineWidth',2);
end
set(gca, 'YLim', [0 10])
xlabel('Electrode Size');
ylabel('Phosephene size (deg)');
set(gca,'FontSize',24);


