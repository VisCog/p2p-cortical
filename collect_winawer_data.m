function collect_winawer_data()
% code hacked from their original code plotting Figure 4 for the following paper:
%  
%   Winawer and Parvizi (2016). Linking Electrical Stimulation of Human
%   Primary Visual Cortex, Size of Affected Cortical Area, Neuronal
%   Responses, and Subjective Experience Neuron. 92(6):1213?1219
%   http://dx.doi.org/10.1016/j.neuron.2016.11.008


sites = 1:5;

opts = getPlotOpts();  % Plot options

D = getData(sites);    % Data per trials

D = deriveData(D, opts); % Some additional derived parameters

% Figure 4a
figure4a(D, opts); 

% Figure 4b
figure4b(D, opts); 

% Figure 4c
figure4c(D, opts); 

% Figure 4d
figure4d(D, opts);

% Figure 4e
figure4e(D, opts); 

end


function opts = getPlotOpts

fieldsToPlot = {'chargedensity' 'chargedensityPerPulse' 'totalCharge' 'chargedensityPerTime'};
opts.fieldToPlot = fieldsToPlot{3};

% Plot cortical area using individual maps ('area') or standard cm function
% ('cm_area')
opts.area = 'area';

% Plot colors
[opts.colors, opts.sites] = getColors;

% Legend text
opts.leg_txt = cellstr(num2str(opts.sites'));

% Axis scale and range
opts.yscale = 'log';
opts.xscale = 'log'; % 'linear'
opts.yl     = 10.^[-1 3]; % 'linear'

switch opts.xscale
    case 'log', opts.xl = 10.^[1.5 4.5];
    case 'linear', opts.xl = [0 1.2e4];
end


end

function fH = figure4b(D, opts)

fH = figure; 
set(fH, 'Color', 'w' ,'name', 'Winawer and Parvizi Figure 4B', 'NumberTitle', 'off');

set(gca, 'FontSize', 30); hold on

x_all = []; for ii = 1:length(D); x_all = [x_all D(ii).stimulation_data.val]; end
x_all = unique(x_all);

for ii = 1:length(D)
    plot(D(ii).stimulation_data.val, D(ii).surface_area.val,...
        'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
end

for ii = 1:length(D)
    
    idx = isfinite(D(ii).surface_area.val);
    x = D(ii).stimulation_data.val(idx)';
    y = D(ii).surface_area.val(idx);
    
    [f, gof] = fit(x,y,'b*x^m', 'StartPoint',[mean(x) / mean(y) 1]);
    pred_y = f(x_all);
    
    disp(f)
    disp(gof)
    
    plot(D(ii).stimulation_data.val, D(ii).surface_area.val, ...
        'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
    plot(x_all,  pred_y, 'Color', opts.colors(ii,:), 'LineWidth', 2)
    
    
end

xlabel('Charge Deposited Per Trial (µC)');
ylabel('Cortical Area (mm^2)');

set(gca, 'YScale', opts.yscale, 'XScale', opts.xscale, 'XLim', 10.^[0 3], 'XTick', 10.^[0 1 2 3])
if strcmp(opts.yscale, 'log'), set(gca, 'YLim', 10.^[-3 3]); end
% if strcmp(xscale, 'log'), set(gca, 'XLim', 10.^[1.5 4.5]); end

xl = get(gca, 'XLim');
plot(xl, 1.13 * [1 1], 'k--')
plot(xl, 4.15 * [1 1], 'k--')

end

function fH = figure4e(D, opts)

x1 = 'chargePerPulse';
x2 = 'frequency'; % 'num_pulses'
% dvs = {'surface_area_indiv' 'surface_area_cm' 'subjective_rating'}; % 'phosphene_area'
dvs = { 'surface_area' 'subjective_rating'};
x1all = []; x2all = [];
for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    x1all = [x1all D(ii).(x1).val(idx)]; 
    x2all = [x2all D(ii).(x2).val(idx)];
end



fH = figure; pos = get(fH, 'Position'); pos(3:4) = [800 600];
set(gcf, 'Color', 'w', 'name', 'Winawer & Parvizi, Figure S4',...
    'NumberTitle', 'off', 'Position', pos);

for z = 1:length(dvs)
    dv = dvs{z};
    
    for ii = 1:length(D)
        
        if all(isnan(D(ii).(dv).val)), skip = true; else, skip = false; end
        
        disp(dv)
        subplot(length(D),length(dvs),(ii-1)*length(dvs)+z)
        set(gca, 'FontSize', 15); hold on;
        if ~skip
            idx = isfinite(D(ii).which_drawing.val);
            
            lm = fitlm([D(ii).(x1).val(idx)', D(ii).(x2).val(idx)'],...
                D(ii).(dv).val(idx)', 'linear',  'varnames', {x1, x2, dv});
            
            lm.plotEffects; xl = get(gca, 'XLim'); set(gca, 'XLim', [-1 1] * xl(2));
        else
            axis off;
        end
        
        if ii == 1, title(dv, 'interpreter', 'none'); end
        
    end
    
    %     subplot(length(D)+1,length(dvs),length(D)*length(dvs)+z)
    %     dvall = []; for ii = 1:5; dvall = [dvall D(ii).(dv)']; end
    %     lm = fitlm([x1all; x2all]', dvall', 'linear');
    %     lm.plotEffects; xl = get(gca, 'XLim'); set(gca, 'XLim', [-1 1] * xl(2));
    
end


fH(2) = figure; pos = get(fH(2), 'Position'); 

set(fH(2), 'Position', [pos(1) pos(2) 400 800], 'Color', 'w', ...
    'name', 'Winawer and Parvizi Figure 4E', 'NumberTitle', 'off');

for z = 1:length(dvs)
    
    xmx = 10.^ceil(log10(max(x1all)));
    xmn = 10.^floor(log10(min(x1all)));
    
    %xl = [8 400];  yl = [4  120];
    xl = [xmn xmx];  yl = [4  120];
    xl(1) = .1;
    xt = 10.^(log10(xmn):log10(xmx)); yt = [10 100];
    %xt = [10 100]; yt = [10 100];
    dv = dvs{z};
    
    for ii = 1:length(D)
        
        idx = isfinite(D(ii).which_drawing.val);
        if all(isnan(idx)), skip = true; else, skip = false; end
        
        subplot(length(D),length(dvs),(ii-1)*length(dvs)+z)
        set(gca, 'FontSize', 20); hold on;
        if ii == 1, title(dv, 'interpreter', 'none'); end
        
        if ~skip
            sz = D(ii).(dv).val(idx);
            sz = sz * 300 / max(sz);
            sz(sz == 0) = eps;
            
            scatter(D(ii).(x1).val(idx), D(ii).(x2).val(idx), sz, ...
                'MarkerFaceColor', opts.colors(ii,:),...
                'MarkerEdgeColor', 'k', 'LineWidth', 1)
            scatter(D(ii).(x1).val(idx), D(ii).(x2).val(idx), sz,...
                'MarkerEdgeColor', 'k', 'LineWidth', 2)
            if ii == length(D)
                xlabel(sprintf('Charge per pulse\n%s', D(ii).chargePerPulse.units));
            end
            if z == 1 && ii == 3, ylabel('Frequency (Hz)'); end
            
            axis([xl yl])
            set(gca, 'YScale', opts.yscale, 'XScale', opts.xscale, 'XTick', xt,'YTick', yt)
            axis square
            plot(xl, [15 15], 'k--', 0.7*[1 1], yl, 'k--')
        end
    end
    
    
end




end

function fH = figure4a(D, opts)
fH = figure; 

set(fH, 'Color', 'w' ,'name', 'Winawer and Parvizi Figure 4A', 'NumberTitle', 'off');

set(gca, 'FontSize', 30); hold on

fit_type = 'power';

x_all = []; 
for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    x_all = [x_all D(ii).stimulation_data.val(idx)]; 
end

x_all = sort(unique(x_all));

for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    switch fit_type
        case 'linear'
            b = regress(D(ii).poly_area.val(idx)', ...
                D(ii).stimulation_data.val(idx)');
            pred_y = x_all * b;
            
            
        case 'power'
            [f, gof] = fit(D(ii).stimulation_data.val(idx)', ...
                D(ii).poly_area.val(idx)','b*x^m', ...
                'StartPoint',[mean(D(ii).poly_area.val(idx)) ...
                / mean(D(ii).stimulation_data.val(idx)) 1]);
            pred_y = f(x_all);
            disp(f)
            disp(gof)
            
    end
    plot(D(ii).stimulation_data.val(idx), D(ii).poly_area.val(idx),...
        'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
    plot(x_all,  pred_y, 'Color', opts.colors(ii,:), 'LineWidth', 2)
end

set(gca, 'YScale', opts.yscale, 'XScale', opts.xscale, 'XTick', 10.^[0 1 2 3], 'XLim', 10.^[0 3])
if strcmp(opts.yscale, 'log'), set(gca, 'YLim', 10.^[-3 3]); end

xlabel('Charge Deposited per Trial (µC)')
ylabel('Phosphene size (deg^2)')



end

function fH = figure4c(D, opts)
fH = figure;
set(gcf, 'Color', 'w' ,'name', 'Winawer and Parvizi Figure 4C');

% all channels, one plot
hold on, set(gca, 'FontSize', 30)
fit_type = 'power';

for ii = 1:numel(D)
    
    [x{ii}, inds] = sort(D(ii).poly_area.val');
    y{ii} = D(ii).surface_area.val(inds);
    idx = isfinite(x{ii});
    x{ii} = x{ii}(idx);
    y{ii} = y{ii}(idx);
    sz = D(ii).stimulation_data.val;
    sz = sz(inds) / max(sz)*200;
    
    
    plot(x{ii}, y{ii}, 'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
    %scatter(x{ii}, y{ii}, sz, 'ko','MarkerFaceColor', opts.colors(ii,:))
    
    switch fit_type            
        case 'power'
            % power law fit
            [f, gof] = fit(x{ii},y{ii},'b*x^m', 'StartPoint', [mean(y{ii})/mean(x{ii}) 1]);
            xpred{ii} = [min(x{ii})/10; x{ii}; max(x{ii})*10];
            ypred{ii} = f(xpred{ii});
            disp(f)
            disp(gof)
        case 'linear'
            % linear fit
            [b{ii}, ~ ,~, ~, stats] = regress(y{ii}, x{ii});
            ypred{ii} = x{ii}* b{ii};
            r2(ii)=(stats(1));
    end
end
for ii = 1:length(D)     
    plot(xpred{ii}, ypred{ii}, '-', 'Color', opts.colors(ii,:), 'LineWidth', 2);
end
xl =  10.^[-3 3];
yl = 10.^[0 3];
set(gca, 'YScale', opts.yscale,  'XScale', 'log', 'YLim', yl, 'XLim',xl, 'XTick', 10.^[-2 0 2])
ylabel('Cortical area (mm^2)'), xlabel('Phosphene size (deg^2)')
legend(opts.leg_txt, 'Location', 'Best')


end

function fH = figure4d(D, opts)

fH = figure; set(fH, 'Color', 'w' ,'name', 'Winawer and Parvizi Figure 4D');

x_all = []; 
for ii = 1:length(D)
    idx = isfinite(D(ii).which_drawing.val);
    x_all = [x_all D(ii).stimulation_data.val(idx)]; 
end

x_all = sort(unique(x_all));


set(gca, 'FontSize', 30); hold on
% title('Subjective intensity rating')

fit_type = 'power';

% % if a trial has the same electrode and same condition number as the trial
% % before, then it must indicate a second drawing on the same trial (subject
% % drew two phosphenes for one stimulation and we only use the first)
% isrepeat =  [1  diff(electrode)] == 0 & [1 diff(condition)] == 0;

for ii = 1:length(D)
    
    if all(isnan(D(ii).subjective_rating.val))
    else
        
        idx = isfinite(D(ii).subjective_rating.val);
        switch fit_type
            case 'linear'
                b = regress(D(ii).subjective_rating.val(idx)', ...
                    D(ii).stimulation_data.val(idx)');
                pred_y = x_all * b;
                
                
            case 'power'
                [f, gof] = fit(D(ii).stimulation_data.val(idx)', ...
                    D(ii).subjective_rating.val(idx)','b*x^m', ...
                    'StartPoint',[mean(D(ii).subjective_rating.val(idx)) ...
                    / mean(D(ii).stimulation_data.val(idx)) 1]);
                pred_y = f(x_all);
                disp(f)
                disp(gof)
        end
        
        plot(D(ii).stimulation_data.val(idx), D(ii).subjective_rating.val(idx),...
            'ko','MarkerFaceColor', opts.colors(ii,:), 'MarkerSize', 12)
        plot(x_all,  pred_y, 'Color', opts.colors(ii,:), 'LineWidth', 2)
        
    end
end

yscale = 'linear';
set(gca, 'YScale', yscale, 'XScale', opts.xscale, 'YLim', [0 11], ...
    'YTick', 0:2:10, 'XLim', 10.^[0 3], 'XTick', 10.^[0 1 2 3])
if strcmp(yscale, 'log'), set(gca, 'YLim', 10.^[0 1]); end

xlabel('Charge Deposited per Trial (µC)')
ylabel('Subjective rating')

end

function D = getData(sites)

pth = fullfile(ebsRootPath, 'data', 'ebs');
for ii = sites
    fname = sprintf('trial_data_site%d', ii);
    D(ii) = load(fullfile(pth, fname));    
end

%[~,~,~,units] = getConditionsFromFile('jt', pth.data);

end

function D = deriveData(D, opts)

for ii = 1:length(D)
    
    % Which stimulation parameter to plot?
    D(ii).stimulation_data.val      = D(ii).(opts.fieldToPlot).val;
    
    % Which measure of surface area to plot? (derived from indivual
    % retinoptic map or from standard CM function)
    D(ii).surface_area.val = D(ii).(opts.area);
    
    % Derive number of pulses from frequency and duration
    D(ii).num_pulses.val            = D(ii).frequency.val .* D(ii).duration.val;
    
    % Derive one rating pre trial from separate ratings for color, motion,
    % brightness
    D(ii).subjective_rating.val     = nanmedian([...
        D(ii).motion.val;...
        D(ii).color.val; ...
        D(ii).brightness.val...
        ]);
    
end

% If we are plotting the surface area using indivual retinotopic maps, then
% we interpolate the area based on the standard CMF. This is because for
% small phosphenes, there may be no voxel whose center is insider the
% phosphene, yet we know the surface area cannot be 0.
for ii = 1:numel(D)
    idx = isfinite(D(ii).which_drawing.val);
    lm = fitlm(D(ii).cm_area.val(idx), D(ii).area.val(idx), 'Intercept', false);    
    D(ii).area.val(idx) = lm.predict(D(ii).cm_area.val(idx)');
    D(ii).surface_area.val = D(ii).(opts.area).val';    
end

end
