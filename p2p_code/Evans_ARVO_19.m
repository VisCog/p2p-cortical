% Evans_ARVO_2019.m
% 
% Generates images used for ARVO 2019 talk

clear v c

v.drawthr = .015;

% Phosphene centers, grabbed from figure 1
xy=[-11.3697    -0.1060
    -13.2493    -3.1039
    -11.1232    -4.7089
    -5.1765     -3.8913
    -5.9160     -4.7089
    -3.5,       -3.7
    0.0         6.2533
    0.0        -2.6497
    -8.9356      2.9828
    0.0         4.4666];

%  xy(:,1) = xy(:,1);
electrodeID = [1,7,13,34,35,41,48,57,60,64];

electrodeList= [1:10];

colorList = hsv(length(electrodeList)+1);

for ii = electrodeList %length(electrodeID)
    v.e(ii).ang =180*atan2(xy(ii,2),xy(ii,1))/pi;
    v.e(ii).ecc = sqrt(xy(ii,1)^2 + xy(ii,2)^2);
    c.e(ii).radius = 2;  % don't know this
end
v.retinaSize = [44,22];v.pixperdeg = 10;
v.retinaCenter = [-10,0];

c.cortexSize = [80,100];
c.cortexCenter = [30,0];


c = p2p_c.define_cortex(c);
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);
c = p2p_c.define_electrodes(c, v);
tp.scFac = 1/10;  %1/10
tp = p2p_c.define_temporalparameters(tp);
c = p2p_c.generate_ef(c);
dx = c.x(1,2)-c.x(1,1);
%%
% show phosphene locations in retinal space
figure(1)
clf
p2p_c.plotretgrid(1, v, gray(64), 1, ['']);
hold on
plot([v.e.ecc].*cos(pi*[v.e.ang]/180),[v.e.ecc].*sin(pi*[v.e.ang]/180),'ko','MarkerFaceColor','w','MarkerSize',14);

for i=electrodeList
    text(xy(i,1),xy(i,2),num2str(electrodeID(i)),'Color','r','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
end
%%
% Generate cortical electrical field composition
cImg = zeros(size(c.e(1).ef));
for ii=electrodeList
    cImg = cImg + c.e(ii).ef;
end
    
cImg = cImg/max(cImg(:));
figure(2)
clf

p2p_c.plotcortgrid(cImg*256, c, gray(256), 2,'');
set(gca,'FontSize',12);
hold on
for i=electrodeList
    tmpz = p2p_c.c2v(c, -xy(i,1)-sqrt(-1)*xy(i,2));    
    text(real(tmpz),imag(tmpz),num2str(electrodeID(i)),'Color','w','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
end

%%
% Generate image of phosphene composition

vImg = zeros(size(v.X));

for ii=electrodeList
    
    %  Draw a phosphene for each electrode
    tmp = [];
    v = p2p_c.generate_rfmap(c, v, ii);
    
    tmp.expname = 'Evans';
    
    tmp.e = ii; % which electrode
    tmp = p2p_c.define_trial(tp,tmp);
    tmp.amp = 5000;
    
    tmp = p2p_c.generate_phosphene(v, tp, tmp);
    vImg = vImg+tmp.maxphos;
end

%%
% Show the phosphene composition image
vImg = vImg/max(vImg(:));

figure(3)
clf
p2p_c.plotretgrid(vImg(:,:,1)*256, v, gray(256), 3,'');
hold on
for i=electrodeList
    text(xy(i,1),xy(i,2),num2str(electrodeID(i)),'Color','b','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
end
%%
% Perceived brightness as a function of stimulation amplitude

ampList = 1000*exp(linspace(log(.5),log(4),7));
electrodeList =[1];
trl ={};
sim_sizes = [];
sim_brightness = [];
for ii=electrodeList
    
    %  Draw a phosphene for each electrode
    tmp = [];
    v = p2p_c.generate_rfmap(c, v, ii);
    
    tmp.expname = 'Evans';
    
    tmp.e = ii; % which electrode
    tmp = p2p_c.define_trial(tp,tmp);
    tmp.amp = 5000;
    
    tmp = p2p_c.generate_phosphene(v, tp, tmp);
    plotcortgrid(c.e(ii).ef * 256, c, gray(256), ii,['subplot(2,1,2); title(''electric field'')']);
%     
    disp(sprintf('Electrode %d of %d',ii,length(v.e)));
    tmp = [];
    tmp.expname = 'Evans';
    for tt=1:length(ampList)
        tmp.amp = ampList(tt);
        tmp.e = ii;
        tmp = p2p_c.define_trial(tp,tmp);
        
        tmp = p2p_c.generate_phosphene(v, tp, tmp);
        tmp.sim_radius= mean([tmp.ellipse(1).sigma_x tmp.ellipse(1).sigma_y]);
        tmp.sim_diameter = 2 * tmp.sim_radius;
    
        tmp.maxresp = max(tmp.resp);
        
        sim_sizes(ii,tt) =  tmp.sim_diameter;
        sim_brightness(ii,tt) =  tmp.sim_brightness;
        
        trl{ii,tt} = tmp;
    end
end




%%
figure(4)

clf
for ii =electrodeList
    %subplot(5,2,ii);
    hold on
    x = log([min(ampList/1000),max(ampList/1000)]);
    p = polyfit(log(ampList/1000),sim_brightness(ii,:),1);
    y = polyval(p,x)
    plot(x,y,'-','Color','k','LineWidth',1);
    plot(log(ampList/1000),(sim_brightness(ii,:)'),'ko','MarkerFaceColor',[.5,.5,.5],'MarkerSize',12);
    
    set(gca,'XTick',log(ampList(1:2:end)/1000));
    set(gca,'YLim',[min(sim_brightness(:)),max(sim_brightness(:))*1.1]);
    set(gca,'YTick',[]);
    xlabel('Current (\muA)')
    ylabel('Brightness');
    text(x(1),max(sim_brightness(:)*1.08),sprintf('E %d',electrodeID(ii)),'FontSize',16);
    logx2raw
    set(gca,'FontSize',16);
    %logy2raw
end

%% Brightnes pulse width curves

pwList = [.15,.25,.5,1,2]/1000;

trl ={};
sim_sizes = [];
sim_pw_brightness = [];
figure(15)
clf
for ii=electrodeList
    
    %  Draw a phosphene for each electrode
    tmp = [];
    v = p2p_c.generate_rfmap(c, v, ii);
    
    tmp.expname = 'Evans';
    
    
    tmp.e = ii; % which electrode
    tmp = p2p_c.define_trial(tp,tmp);
    tmp.amp = 5000;
    
    tmp = p2p_c.generate_phosphene(v, tp, tmp);
    
    disp(sprintf('Electrode %d of %d',ii,length(v.e)));
    tmp = [];
    tmp.expname = 'Evans';
    for tt=1:length(pwList)
        tmp.e = ii;
        tmp.amp = 3000;
        tmp.pw = pwList(tt);
        tmp = p2p_c.define_trial(tp,tmp);
        
        tmp = p2p_c.generate_phosphene(v, tp, tmp);
        tmp.sim_radius= mean([tmp.ellipse(1).sigma_x tmp.ellipse(1).sigma_y]);
        tmp.sim_diameter = 2 * tmp.sim_radius;
    
        tmp.maxresp = max(tmp.resp);
        
        sim_sizes(ii,tt) =  tmp.sim_diameter;
        sim_pw_brightness(ii,tt) =  tmp.sim_brightness;
        
        trl{ii,tt} = tmp;
        figure(15)
        hold on
        plot(tmp.resp);
    end
end

%%
figure(5)

clf
for ii =electrodeList
%     subplot(5,2,ii);
%     hold on
    x = linspace(log(min(pwList*1000)),log(max(pwList*1000)),21);
    p = polyfit(log(pwList*1000),sim_pw_brightness(ii,:),2);
    y = polyval(p,x);
    plot(x,y,'-','Color',colorList(ii,:),'LineWidth',2);
    plot(log(pwList*1000),(sim_pw_brightness(ii,:)'),'ko','MarkerFaceColor',colorList(ii,:));
    
    set(gca,'XTick',log(pwList*1000));
    set(gca,'YLim',[min(sim_pw_brightness(:)),max(sim_pw_brightness(:))*1.1]);
    set(gca,'YTick',[]);
    xlabel('Pulse Width (ms)')
    ylabel('Brightness');
    text(x(1),max(sim_pw_brightness(:)),sprintf('E %d',electrodeID(ii)));
    logx2raw
    %logy2raw
end
return

%% Brightness frequency curves

freqList = [12.5,50,100,200,400];
trl ={};
sim_sizes = [];
sim_freq_brightness = [];
figure(20)
clf
for ii=electrodeList
    
    %  Draw a phosphene for each electrode
    tmp = [];
    v = p2p_c.generate_rfmap(c, v, ii);
    
    tmp.expname = 'Evans';
    
    
    tmp.e = ii; % which electrode
    tmp = p2p_c.define_trial(tp,tmp);
    tmp.amp = 5000;
    
    tmp = p2p_c.generate_phosphene(v, tp, tmp);
    
    disp(sprintf('Electrode %d of %d',ii,length(v.e)));
    tmp = [];
    tmp.expname = 'Evans';
    for tt=1:length(freqList)
        tmp.e = ii;
        tmp.amp = 5000;
        %  tmp.pw = 25/1000;
        tmp.dur = .5;
        tmp.freq = freqList(tt);
        tmp = p2p_c.define_trial(tp,tmp);
        
        figure(20)
        hold on
        plot(tmp.t,tmp.pt(1:length(tmp.t)));
        
        tmp = p2p_c.generate_phosphene(v, tp, tmp);
        tmp.sim_radius= mean([tmp.ellipse(1).sigma_x tmp.ellipse(1).sigma_y]);
        tmp.sim_diameter = 2 * tmp.sim_radius;
        

        
        
        tmp.maxresp = max(tmp.resp);
        sim_sizes(ii,tt) =  tmp.sim_diameter;
        sim_freq_brightness(ii,tt) =  tmp.sim_brightness;
        
        trl{ii,tt} = tmp;
        
    end
end

%%
figure(6)

sim_freq_brightness = log(sim_freq_brightness);
clf
for ii =electrodeList
%     subplot(5,2,ii);
%     hold on
    x = linspace(log(min(freqList)),log(max(freqList)),21);
    p = polyfit(log(freqList),sim_freq_brightness(ii,:),2);
    y = polyval(p,x);
    plot(x,y,'-','Color',colorList(ii,:),'LineWidth',2);
    plot(log(freqList),(sim_freq_brightness(ii,:)'),'ko','MarkerFaceColor',colorList(ii,:));
    
    set(gca,'XTick',log([12.5,25,50,100,200,400]));
    set(gca,'YLim',[min(sim_freq_brightness(:)),max(sim_freq_brightness(:))*1.1]);
    set(gca,'YTick',[]);
    xlabel('Frequency (Hz)')
    ylabel('Brightness');
    text(x(1),max(sim_freq_brightness(:)),sprintf('E %d',electrodeID(ii)));
    logx2raw
    %logy2raw
end

