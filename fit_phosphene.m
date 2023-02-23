function [err, p] = fit_phosphene(p,trl)
tmp = mean(trl.maxphos, 3);
phos = tmp(round(size(tmp,1)/2), :);
pred_phos =max(phos).*normpdf(1:length(phos), p.mu, p.sigma );
err = sum((pred_phos-phos).^2);
end