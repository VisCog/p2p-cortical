function err = findDrawingThreshold(p, sim_phosphene, respfac, real_area)

% finds the threshold that minimizes difference between simulated and real area drawings across all conditions
for cc = 1: length(respfac)
    sim_area(1,cc) = length( find(sim_phosphene * respfac(cc)>p.thr) );
end

err = sum(sqrt(((log10(sim_area)-log10(real_area)).^2)));
end
