function err = findDrawingThreshold(p, sim_phosphene, respfac, real_area)

% finds the threshold that minimizes difference between simulated and real area drawings across all conditions
for cc = 1: length(respfac)
    sim_area(1,cc) = sum((sim_phosphene * respfac(cc)>p.thr));
end

err = sum(sqrt(((sim_area-real_area).^2)));
end
