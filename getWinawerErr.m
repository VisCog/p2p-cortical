function [err,S]=  getWinawerErr(p,v_electrode,c_electrode,tp,site,W)

colorList =  [ 0.5 0 0 ;0.5 1 0.5; 1  0.8125 0;0 0.8750  1; 0  0 1.0000];

tp.scFac = p.scFac;
err =0;
clear S
for ii=1:length(site)
    v = v_electrode(ii);
    c = c_electrode(ii);
    
    v.drawthr = p.drawthr;
    clear trl
    goodVals = ~isnan(W{ii}.area) & W{ii}.area>0;
    S{ii}.area = zeros(size(W{ii}.area));
    for tt = find(goodVals)
        clear tmp
        tmp = site(ii).trl(tt);
        tmp.e = ii; % which electrode
        tmp = p2p_c.define_trial(tp,tmp);
        
        tmp = p2p_c.generate_phosphene(v, tp, tmp);
        
        S{ii}.area(tt) = tmp.sim_area;
        %         
        %         err = err + sum((log(S{ii}.area(goodVals)+da)-log(W{ii}.area(goodVals)+da)).^2);
    end
   % err = err + sum((sqrt(S{ii}.area(goodVals))-sqrt(W{ii}.area(goodVals))).^2);
   da = .000001;
    err = err + sum((log(S{ii}.area(goodVals)+da)-log(W{ii}.area(goodVals)+da)).^2);
end
figure(1)
clf
subplot(1,2,1)
hold on
for ii=1:length(site)
    plot(log(W{ii}.charge),log(W{ii}.area),'ko','MarkerFaceColor',colorList(ii,:));

    set(gca,'XTick',log([1,10,100,1000]));
    set(gca,'XLim',[log(1),log(1000)]);
    set(gca,'YTick',log([.01,.1,1,10,100,1000]));
    set(gca,'YLim',[log(.01),log(1000)]);
    
    
    logx2raw(exp(1),0);
    logy2raw(exp(1),0);
    xlabel('Charge Deposited per Trial (\muC)');
    ylabel('Phosphene area (deg^2)');
end 

subplot(1,2,2)
hold on
for ii=1:length(site)
    plot(log(W{ii}.charge),log(S{ii}.area),'ko','MarkerFaceColor',colorList(ii,:));

    set(gca,'XTick',log([1,10,100,1000]));
    set(gca,'XLim',[log(1),log(1000)]);
    set(gca,'YTick',log([.01,.1,1,10,100,1000]));
    set(gca,'YLim',[log(.01),log(1000)]);
    
    
    logx2raw(exp(1),0);
    logy2raw(exp(1),0);
    xlabel('Charge Deposited per Trial (\muC)');
    ylabel('Phosphene area (deg^2)');
end 
drawnow
figure(2)
clf
hold on
for ii=1:length(site)
    plot(log(W{ii}.area),log(S{ii}.area),'ko','MarkerFaceColor',colorList(ii,:));

    set(gca,'XTick',log([.01,.1,1,10,100,1000]));
    set(gca,'XLim',[log(.01),log(1000)]);
    set(gca,'YTick',log([.01,.1,1,10,100,1000]));
    set(gca,'YLim',[log(.01),log(1000)]);
    
    
    logx2raw(exp(1),0);
    logy2raw(exp(1),0);
    xlabel('Real Phosphene area (deg^2)');
    ylabel('Simulated Phosphene area (deg^2)');
    title(sprintf('drawthr = %5.4f, scFac = %5.4f, err = %5.4f',p.drawthr,p.scFac,err));
    grid
    axis equal
end 
drawnow



disp(sprintf('drawthr = %5.4f, scFac = %5.4f, err = %5.4f',p.drawthr,p.scFac,err));



