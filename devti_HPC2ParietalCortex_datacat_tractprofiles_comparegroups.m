% Read in tractprofiles and plot tractprofiles for each age group
% and then plot mean and standard deviations for age group.

clear all; close all; clc
format shortG

blprojectid = 'proj-5e9ef616f1745d93d3f6b992';

wm = {'fa', 'md'};

fontname = 'Arial';
fontsizex = 16; fontsizey = 12;
fontangle = 'italic';
fontsmoothing = 'off';
yticklength = 0;
xticklength = 0.05;
linewidth = .15;
alpha = .3;
save_figures = 'yes';

% Set working directories.
rootDir = '/Volumes/240/devti_HPC2ParietalCortex/';
addpath(genpath(fullfile(rootDir, 'proj-5e9ef616f1745d93d3f6b992/')));

% Get age_group.
load(fullfile(rootDir, 'supportFiles/data.mat'))
bdata = array2table(data, 'VariableNames', {'subID', 'cov_age', 'iq', 'gp_age', 'a', 'b', 'c', 'd', 'e'});

% Add fix so to account for leading zeros in file names.
% for k = 1:length(bdata.subID)
%
%     if bdata.subID(k) < 10
%
%         subID{k, 1} = strcat('00', num2str(bdata.subID(k)));
%
%     else
%
%         subID{k, 1} = strcat('0', num2str(bdata.subID(k)));
%
%     end
%
% end

% Identify outliers to be removed.
outlier = []; %[128 315 318]; %NOTE: need to update to not plot the outliers

fcount = 0;
for w = 2%:length(wm)
    
    wm_measure = wm{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    % Load in each tract's tractography measures for this subject.
    sub_count = 0;
    for i = 1:size(grp_contents, 1)
        
        % Display current sub ID.
        disp(grp_contents(i).name)
        
        % Update subject counter for when not all subjects are used/needed.
        sub_count = sub_count + 1;
        
        %% LEFT HEMISPHERE
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, 'dt-neuro-tractprofile*', 'profiles', '*left*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
        
        if ~isempty(sub_contents_tractprofiles)
            
            % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
            if i == 1
                
                tract_LH = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
                
            end
            
            % Read in data for this subject and this tract.
            data_temp = readtable(fullfile(sub_contents_tractprofiles(1).folder, sub_contents_tractprofiles(1).name));
            
            % Get middle 80%.
            start = size(data_temp, 1)*.1;
            stop = size(data_temp, 1)*.9;
            
            % Read in mean WM measure.
            if strcmp(wm_measure, 'fa')
                
                m_LH(:, sub_count) = data_temp.fa_1(start:stop);
                sd_LH(:, sub_count) = data_temp.fa_2(start:stop);
                
            elseif strcmp(wm_measure, 'md')
                
                m_LH(:, sub_count) = data_temp.md_1(start:stop);
                sd_LH(:, sub_count) = data_temp.md_2(start:stop);
                
            end
            
            % Grab tract name for grouping variable.
            tract_LH(:, sub_count) = {sub_contents_tractprofiles(1).name(1:end-13)};
            
            % Grab subID.
            sub_LH(sub_count) = str2num(grp_contents(i).name(end-2:end));
            
            % Grab age.
            age_LH(sub_count) = bdata.cov_age(bdata.subID == sub_LH(sub_count));
            
            % Grab age group.
            gp_LH(sub_count) = bdata.gp_age(bdata.subID == sub_LH(sub_count));
            
            clear data_temp sub_contents_tractprofiles
            
        end % if
        
        %     end % end sub_contents
        
        %% RIGHT HEMISPHERE
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, 'dt-neuro-tractprofile*', 'profiles', '*right*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
        
        % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
        if i == 1
            
            tract_RH = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
            
        end
        
        % Read in data for this subject and this tract.
        data_temp = readtable(fullfile(sub_contents_tractprofiles(1).folder, sub_contents_tractprofiles(1).name));
        
        % Get middle 80%.
        start = size(data_temp, 1)*.1;
        stop = size(data_temp, 1)*.9;
        
        % Read in mean WM measure.
        if strcmp(wm_measure, 'fa')
            
            m_RH(:, sub_count) = data_temp.fa_1(start:stop);
            sd_RH(:, sub_count) = data_temp.fa_2(start:stop);
            
        elseif strcmp(wm_measure, 'md')
            
            m_RH(:, sub_count) = data_temp.md_1(start:stop);
            sd_RH(:, sub_count) = data_temp.md_2(start:stop);
            
        end
        
        % Grab tract name for grouping variable.
        tract_RH(sub_count) = {sub_contents_tractprofiles(1).name(1:end-13)};
        
        % Grab subID.
        sub_RH(sub_count) = str2num(grp_contents(i).name(end-2:end));
        
        % Grab age group.
        age_RH(sub_count) = bdata.cov_age(bdata.subID == sub_RH(sub_count));
        
        % Grab age group.
        gp_RH(sub_count) = bdata.gp_age(bdata.subID == sub_RH(sub_count));
        
        clear data_temp sub_contents_tractprofiles
        
    end % group_contents
    
    %% PLOT LEFT HEMISPHERE
    
    % Find empty cells and fill with 'empty'.
    t = find(cellfun(@isempty,tract_LH));
    tract_LH(t) = {'empty'};
    
    % Get a list of unique tract names.
    list_tract = unique(tract_LH);
    
    % Get a list of unique sub IDs.
    subID = unique(sub_LH(sub_LH~=0));
    
    if ~strcmp(list_tract{1}, 'empty')
        
        % Open a new figure and hold.
        figure(w)
        hold on;
        
        % Set up legend.
        scatter([-9 -10], [-9 -10], 'Marker', 'o', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerEdgeColor', [0.6350 0.0780 0.1840])
        scatter([-9 -10], [-9 -10], 'Marker', 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', [0 0.4470 0.7410])
        scatter([-9 -10], [-9 -10], 'Marker', 'o', 'MarkerFaceColor', [0.41176 0.41176 0.41176], 'MarkerEdgeColor', [0.41176 0.41176 0.41176])
        
        count = 0; oc_count = 0; a_count = 0;
        for s = 1:length(subID)
            
            if subID(s) ~= 0 && ~ismember(subID(s), outlier)
                
                % Find entries that are for this subject.
                s_idx = sub_LH == subID(s);
                
                % Subset the thing so that we only plot for this subject.
                t_temp = m_LH(:, find(s_idx));
                
                if ~isempty(t_temp)
                    
                    count = count + 1;
                    
                    % Get age group to colorcode.
                    if gp_RH(count) == 1
                        cmap = [0.6350 0.0780 0.1840]; %red, child
                    elseif gp_RH(count) ==2
                        cmap = [0 0.4470 0.7410]; %blue, adolescent
                    elseif gp_RH(count) == 3
                        cmap = [0.41176 0.41176 0.41176]; %gray, adult
                    end
                    
                    % Plot. Flip to put the HPC first and parietal cortex second.
                    plot(flip(t_temp), 'LineStyle', '-', 'LineWidth', linewidth, 'Color', [cmap alpha])
                    
                    % Collect. Flip to put the HPC first and parietal cortex second.
                    y(:, count) = flip(t_temp);
                    
                end
                
            end % if subID{s} ~= 0
            
        end %sub
        
    end %if strcmp
    
    if strcmp(wm{w}, 'fa')
        
        ylim_lo = 0.10; ylim_hi = 0.70;
        
    elseif strcmp(wm{w}, 'md')
        
        ylim_lo = 0.5; ylim_hi = 1.5;
        
    end
    
    % Plot means and standard deviations: children.
    cmap = [0.6350 0.0780 0.1840]; %red, child
    y_child = y(:, gp_LH == 1);
    plot(nanmean(y_child, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap)
    hi = nanmean(y_child, 2) + nanstd(y_child, 0, 2); lo = nanmean(y_child, 2) - nanstd(y_child, 0, 2);
    x = (1:size(nanmean(y_child, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap);
    set(hp1, 'facecolor', cmap, 'edgecolor', 'none', 'facealpha', .2);
    
    % Plot means and standard deviations: adolescent.
    cmap = [0 0.4470 0.7410]; %blue, adolescent
    y_adolescent = y(:, gp_LH == 2);
    plot(nanmean(y_adolescent, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap)
    hi = nanmean(y_adolescent, 2) + nanstd(y_adolescent, 0, 2); lo = nanmean(y_adolescent, 2) - nanstd(y_adolescent, 0, 2);
    x = (1:size(nanmean(y_adolescent, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap);
    set(hp1, 'facecolor', cmap, 'edgecolor', 'none', 'facealpha', .2);
    
    % Plot means and standard deviations: adult.
    cmap = [0.41176 0.41176 0.41176]; %gray, adult
    y_adult = y(:, gp_LH == 3);
    plot(nanmean(y_adult, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap)
    hi = nanmean(y_adult, 2) + nanstd(y_adult, 0, 2); lo = nanmean(y_adult, 2) - nanstd(y_adult, 0, 2);
    x = (1:size(nanmean(y_adult, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap);
    set(hp1, 'facecolor', cmap, 'edgecolor', 'none', 'facealpha', .2);
    
    % yaxis - left
    ylabel(wm_measure);
    %     ylim_lo = 0.30; ylim_hi = 0.70;
    yax = get(gca,'yaxis');
    yax.Limits = [ylim_lo ylim_hi];
    yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
    yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    %     yax.Direction = 'normal';
    yax.FontName = fontname;
    yax.FontSize = fontsizey;
    yax.FontSmoothing = fontsmoothing;
    
    % xaxis
    xlabel('Location along tract');
    xlim_lo = 0; xlim_hi = 160;
    xax = get(gca, 'xaxis');
    xax.Limits = [xlim_lo xlim_hi];
    xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
    xax.TickDirection = 'out';
    xax.TickLength = [xticklength xticklength];
    xax.TickLabels = {num2str(xlim_lo, '%1.0f'), '', num2str(xlim_hi, '%1.0f')};
    xax.FontName = fontname;
    xax.FontSize = fontsizex;
    xax.FontAngle = fontangle;
    xax.FontSmoothing = fontsmoothing;
    legend({'Children', 'Adolescent', 'Adult'})
    legend box off
    title('Left HPC to Left Parietal Cortex');
    box off;
    
    print(fullfile(rootDir, 'plots', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_LH']), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_LH']), '-depsc')
    
    hold off;
    
    %% PLOT RIGHT HEMISPHERE
    
    % Find empty cells and fill with 'empty'.
    t = find(cellfun(@isempty,tract_RH));
    tract_RH(t) = {'empty'};
    
    % Get a list of unique tract names.
    list_tract = unique(tract_RH);
    
    % Get a list of unique sub IDs.
    subID = unique(sub_RH(sub_RH~=0));
    
    if ~strcmp(list_tract{1}, 'empty')
        
        % Open a new figure and hold.
        figure(w+1)
        hold on;
        
        % Set up legend.
        scatter([-9 -10], [-9 -10], 'Marker', 'o', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerEdgeColor', [0.6350 0.0780 0.1840])
        scatter([-9 -10], [-9 -10], 'Marker', 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', [0 0.4470 0.7410])
        scatter([-9 -10], [-9 -10], 'Marker', 'o', 'MarkerFaceColor', [0.41176 0.41176 0.41176], 'MarkerEdgeColor', [0.41176 0.41176 0.41176])
        
        count = 0;
        for s = 1:length(subID)
            
            if subID(s) ~= 0 && ~ismember(subID(s), outlier)
                
                % Find entries that are for this subject.
                s_idx = sub_RH == subID(s);
                
                % Subset the thing so that we only plot for this subject.
                t_temp = m_RH(:, find(s_idx));
                
                if ~isempty(t_temp)
                    
                    count = count + 1;
                    
                    % Get age group to colorcode.
                    if gp_RH(count) == 1
                        cmap = [0.6350 0.0780 0.1840]; %red, child
                    elseif gp_RH(count) ==2
                        cmap = [0 0.4470 0.7410]; %blue, adolescent
                    elseif gp_RH(count) == 3
                        cmap = [0.41176 0.41176 0.41176]; %gray, adult
                    end
                    
                    % Plot. Flip to put the HPC first and parietal cortex second.
                    plot(flip(t_temp), 'LineStyle', '-', 'LineWidth', linewidth, 'Color', [cmap alpha])
                    
                    % Collect. Flip to put the HPC first and parietal cortex second.
                    y(:, count) = flip(t_temp);
                    
                end
                
            end % if subID{s} ~= 0
            
        end %sub
        
    end %if strcmp
    
    % Plot means and standard deviations: children.
    cmap = [0.6350 0.0780 0.1840]; %red, child
    y_child = y(:, gp_RH == 1);
    plot(nanmean(y_child, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap)
    hi = nanmean(y_child, 2) + nanstd(y_child, 0, 2); lo = nanmean(y_child, 2) - nanstd(y_child, 0, 2);
    x = (1:size(nanmean(y_child, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap);
    set(hp1, 'facecolor', cmap, 'edgecolor', 'none', 'facealpha', .2);
    
    % Plot means and standard deviations: adolescent.
    cmap = [0 0.4470 0.7410]; %blue, adolescent
    y_adolescent = y(:, gp_RH == 2);
    plot(nanmean(y_adolescent, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap)
    hi = nanmean(y_adolescent, 2) + nanstd(y_adolescent, 0, 2); lo = nanmean(y_adolescent, 2) - nanstd(y_adolescent, 0, 2);
    x = (1:size(nanmean(y_adolescent, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap);
    set(hp1, 'facecolor', cmap, 'edgecolor', 'none', 'facealpha', .2);
    
    % Plot means and standard deviations: adult.
    cmap = [0.41176 0.41176 0.41176]; %gray, adult
    y_adult = y(:, gp_RH == 3);
    plot(nanmean(y_adult, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap)
    hi = nanmean(y_adult, 2) + nanstd(y_adult, 0, 2); lo = nanmean(y_adult, 2) - nanstd(y_adult, 0, 2);
    x = (1:size(nanmean(y_adult, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap);
    set(hp1, 'facecolor', cmap, 'edgecolor', 'none', 'facealpha', .2);
    
    % yaxis - left
    ylabel(wm_measure);
    %     ylim_lo = 0.30; ylim_hi = 0.70;
    yax = get(gca,'yaxis');
    yax.Limits = [ylim_lo ylim_hi];
    yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
    yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    %     yax.Direction = 'normal';
    yax.FontName = fontname;
    yax.FontSize = fontsizey;
    yax.FontSmoothing = fontsmoothing;
    
    % xaxis
    xlabel('Location along tract');
    xlim_lo = 0; xlim_hi = 160;
    xax = get(gca, 'xaxis');
    xax.Limits = [xlim_lo xlim_hi];
    xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
    xax.TickDirection = 'out';
    xax.TickLength = [xticklength xticklength];
    xax.TickLabels = {num2str(xlim_lo, '%1.0f'), '', num2str(xlim_hi, '%1.0f')};
    xax.FontName = fontname;
    xax.FontSize = fontsizex;
    xax.FontAngle = fontangle;
    xax.FontSmoothing = fontsmoothing;
    legend({'Children', 'Adolescent', 'Adult'})
    legend box off
    title('Right HPC to Right Parietal Cortex');
    box off;
    
    print(fullfile(rootDir, 'plots', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_RH']), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_RH']), '-depsc')
    
    hold off;
    
    % Get group means and confidence intervals.
    % Use 2*standard error as a close approximate of 95% confidence interval.
    m_child_LH = nanmean(m_LH(gp_LH == 1, 2)); se2_child_LH = 2*(nanstd(m_LH(gp_LH == 1, 2))/sqrt(length(m_LH(gp_LH == 1, 2))));
    m_adolescent_LH = nanmean(m_LH(gp_LH == 2, 2)); se2_adolescent_LH = 2*(nanstd(m_LH(gp_LH == 2))/sqrt(length(m_LH(gp_LH == 2, 2))));
    m_adult_LH = nanmean(m_LH(gp_LH == 3, 2)); se2_adult_LH = 2*(nanstd(m_LH(gp_LH == 3, 2))/sqrt(length(m_LH(gp_LH == 3, 2))));
    
    m_child_RH = nanmean(m_RH(gp_RH == 1, 2)); se2_child_RH = 2*(nanstd(m_RH(gp_RH == 1, 2))/sqrt(length(m_RH(gp_RH == 1, 2))));
    m_adolescent_RH = nanmean(m_RH(gp_RH == 2, 2)); se2_adolescent_RH = 2*(nanstd(m_RH(gp_RH == 2))/sqrt(length(m_RH(gp_RH == 2, 2))));
    m_adult_RH = nanmean(m_RH(gp_RH == 3, 2)); se2_adult_RH = 2*(nanstd(m_RH(gp_RH == 3, 2))/sqrt(length(m_RH(gp_RH == 3, 2))));
    
    capsize = 0;
    marker = 'o';
    cmap_child = [0.6350 0.0780 0.1840]; %red, child
    cmap_adolescent = [0 0.4470 0.7410]; %blue, adolescent
    cmap_adult = [0.41176 0.41176 0.41176]; %gray, adult
    linewidth = 1.5;
    linestyle = 'none';
    markersize = 10;
    fontname = 'Arial';
    fontsize = 16;
    fontangle = 'italic';
    yticklength = 0.05;
    xticklength = 0;
    xtickvalues = [1 2 3 4 5 6];
    xlim_lo = 0.5; xlim_hi = 6.5;
    save_figures = 'yes';
    
    % Plot means. %NOTE: Change errorbar to plot so that I can do whiskers instead of T-bars.
    figure(3)
    errorbars = errorbar(xtickvalues(1), m_child_LH, se2_child_LH, se2_child_LH);
    set(errorbars, 'Color', cmap_child, 'Marker', marker, 'MarkerEdgeColor', cmap_child, 'MarkerFaceColor', cmap_child, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    hold on;
    errorbars = errorbar(xtickvalues(2), m_adolescent_LH, se2_adolescent_LH, se2_adolescent_LH);
    set(errorbars, 'Color', cmap_adolescent, 'Marker', marker, 'MarkerEdgeColor', cmap_adolescent, 'MarkerFaceColor', cmap_adolescent, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    errorbars = errorbar(xtickvalues(3), m_adult_LH, se2_adult_LH, se2_adult_LH);
    set(errorbars, 'Color', cmap_adult, 'Marker', marker, 'MarkerEdgeColor', cmap_adult, 'MarkerFaceColor', cmap_adult, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    errorbars = errorbar(xtickvalues(4), m_child_RH, se2_child_RH, se2_child_RH);
    set(errorbars, 'Color', cmap_child, 'Marker', marker, 'MarkerEdgeColor', cmap_child, 'MarkerFaceColor', cmap_child, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    errorbars = errorbar(xtickvalues(5), m_adolescent_RH, se2_adolescent_RH, se2_adolescent_RH);
    set(errorbars, 'Color', cmap_adolescent, 'Marker', marker, 'MarkerEdgeColor', cmap_adolescent, 'MarkerFaceColor', cmap_adolescent, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    errorbars = errorbar(xtickvalues(6), m_adult_RH, se2_adult_RH, se2_adult_RH);
    set(errorbars, 'Color', cmap_adult, 'Marker', marker, 'MarkerEdgeColor', cmap_adult, 'MarkerFaceColor', cmap_adult, 'MarkerSize', markersize, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    
    % xaxis
    xax = get(gca, 'xaxis');
    xax.Limits = [xlim_lo xlim_hi];
    xax.TickValues = xtickvalues;
    xax.TickDirection = 'out';
    xax.TickLength = [xticklength xticklength];
    xlabels = {'', 'Left Hemisphere', '', '', 'Right Hemisphere', '', ''};
    xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
    xax.TickLabels = xlabels;
    xax.FontName = fontname;
    xax.FontSize = fontsize;
    xax.FontAngle = fontangle;
    
    % yaxis
    yax = get(gca,'yaxis');
    yax.Limits = [ylim_lo ylim_hi];
    yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
    yax.TickDirection = 'out';
    yax.TickLength = [yticklength yticklength];
    yax.TickLabels = {num2str(ylim_lo, '%1.2f'), '', num2str(ylim_hi, '%1.2f')};
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    
    plot([3.5 3.5], [ylim_lo ylim_hi], 'k:')
    
    legend({'Children', 'Adolescents', 'Adults'})
    legend box off
    
    % general
    a = gca;
    %     a.TitleFontWeight = 'normal';
    box off
    
    if strcmp(wm{w}, 'fa')
        a.YLabel.String = 'Fractional Anisotropy (FA)';
    elseif strcmp(wm{w}, 'md')
        a.YLabel.String = 'Mean Diffusivity (MD)';
    end
    a.YLabel.FontSize = fontsize;
    pbaspect([1 1 1])
    
    %     pos=get(gca,'Position');
    %     pos1=pos-[0 .02 0 0];
    %     set(gca,'Position', pos1);
    
    % Write.
    if strcmp(save_figures, 'yes')
        
        print(fullfile(rootDir, 'plots', ['plot_tractprofile_avg_' wm{w} ]), '-dpng')
        %         print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofile_avg_' wm{w} ]), '-depsc')
        
    end
    
    hold off;
    
end % w







