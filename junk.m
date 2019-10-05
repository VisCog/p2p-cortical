figure(10); clf
asymptoteS = 1600;
e50S = 300;% electrical semisaturation constant

asymptoteB = 1500;
slopeB = 0.005;
e50B = 0;

R2 =0:2000;
R3_S = asymptoteS .* (R2.^2./(R2.^2 + e50S.^2));

R3_B =  asymptoteB .*(1-(exp(-slopeB.*(R2-e50B))));

asymptoteC = 1500;
meanC = 750;
sigmoidC = 175;
R3_C =  asymptoteC .* normcdf(R2, meanC, sigmoidC);

plot(R2, R3_S, 'r', R2, R3_B, 'b', R2, R3_C, 'g')
