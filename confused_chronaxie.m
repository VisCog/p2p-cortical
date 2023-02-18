tauList = [.1 1];        % time constant (msec)
ampList = [ 1 2];
     % threshold parameter
cList = ['r', 'b']
d = [0.01:.01:1]; % pulse durations
    figure(1)
    clf
for t = 1:length(tauList)
    tau = tauList(t);
    b = ampList(t);
    p.tau = tau; p.amp = b;
    A = b./(tau*(1-exp(-d/tau)));
    round(min(A))
    round(max(A))
     A2  = p2p_c.chronaxie(p, d);
        plot(d,A2, 'Color', cList(t), 'LineStyle', '--'); hold on
        plot(d,A, 'Color', cList(t)); hold on

end
%set(gca,'YLim',[0,max(A)*1.1]);

xlabel('Pulse Duration (msec)');
ylabel('Threshold mA');