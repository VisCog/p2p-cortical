% p2p_init
% 
% sets up initial model parameters and pulse trains
% then calls a variety of models and compares their output
% 
% written IF 11/16/2017 based on a long, long history of code

%% Model Parameters
% core parameters
p.tau1 = .23/1000;  %.24 - .65
p.tau2_ca = 45.25/1000;  %38-57
p.tau3 =  26.25/1000; % 24-33ex
p.e = 8.73; %2.25;  %2-3 for threshold or 8-10 for suprathreshold

% leak out of charge accumulation
p.flag_cl=0; % 1 if you want to charge to leak back out of the system 
p.tau2_cl = p.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca

% nonlinearity parameters
p.slope=.3; % larger number = shallower slope
p.asymptote=14;
p.shift=47; % shifts curve along x-axis

% models have various flags
%   p.FLAG_cl decides whether or not the charge accumulation leaks away
%   again, 0 no return to baseline, 1 return to baseline. Model was
%   originally developed using 0, but 1 probably correct

% p.FLAG_causal decides whether one scales using a non-causal (as in the
% original model or a causal model (probably correct). 0 = non-causal,
% 1 = causal
% 
%thresholdR2=13.1; % R2 at threshold, need the nonlinearity to be near zero here
%threshold1_5R2=16.37; %
%scFac=p.asymptote./(1+exp(-(p.maxR3-p.shift)./p.slope));

clear STIM


%% Stimulation parameters

STIM.tsample = .01/1000
% Generate tsform vector
STIM.t = 0:STIM.tsample:STIM.dur-STIM.tsample;

dt = STIM.t(2)-STIM.t(1);
n = length(freq);

for i=1   
    sawtooth = freq(i)*mod(STIM.t,1/freq(i));
    on  = sawtooth > STIM.pulsedur*freq(i) & sawtooth < 2*STIM.pulsedur*freq(i);
    off = sawtooth < STIM.pulsedur*freq(i);
    tsform(i,:) = STIM.amp.*(on-off);
end

returnlist={'R4'};
%% original Nanduri model, convolution style
for i=1
    STIM.freq = freq(i);
    STIM.tsform=tsform(i,:);

    p.FLAG_cl=0; p.FLAG_causal=0; % charge only accumulated, non-causal model
    p.maxR3=142.3001;
    out_conv=p2p_devyani_conv(p, STIM, returnlist); % convolution method
    p.maxR3=max(out_conv.R3(:));
    out_fe=p2p_devyani_finite_element(p, STIM, returnlist); % fe method
    
    
    p.FLAG_cl=0;p.FLAG_causal=1; % charge only accumulated, causal model
    out_fe_causal=p2p_devyani_finite_element(p, STIM, returnlist);
    
    p.FLAG_cl=1;p.FLAG_causal=0; % charge accumulates then leaks, non-causal model
    out_conv_cl=p2p_devyani_conv(p, STIM, returnlist);p.maxR3=max(out_conv.R3(:));
    out_fe_cl=p2p_devyani_finite_element(p, STIM, returnlist);
    p.FLAG_cl=1;p.FLAG_causal=1; % charge accumulates then leaks, causal model
    out_fe_cl_causal=p2p_devyani_finite_element(p, STIM, returnlist);
    
end

%% plot data
for r=1:length(returnlist)
    figure(r); clf
    set(gcf, 'Name', returnlist{r});
    h(1)=plot(STIM.t, eval(['out_conv.', returnlist{r}]), 'k-'); hold on
    h(2)=plot(STIM.t, eval(['out_conv_cl.', returnlist{r}]), 'g-'); hold on
    
    h(3)=plot(STIM.t, eval(['out_fe.', returnlist{r}]), 'm--'); hold on
    h(4)=plot(STIM.t, eval(['out_fe_cl.', returnlist{r}]), 'b--'); hold on
    
    h(5)=plot(STIM.t, eval(['out_fe_causal.', returnlist{r}]), 'r:'); hold on
    h(6)=plot(STIM.t, eval(['out_fe_cl_causal.', returnlist{r}]), 'b-.'); hold on
    
    legend(h(1:length(h)),{'out_conv.' 'out_conv_cl' 'out_fe' 'out_fe_cl' 'out_fe_causal' 'out_fe_cl_causal'},  'Interpreter', 'none')
   
    xlabel('time')
    ylabel('response');
end

