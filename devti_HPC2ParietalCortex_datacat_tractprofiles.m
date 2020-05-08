% Read in tractprofiles and plot mean and standard deviations tract.

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
alpha = .5;
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
for w = 1:length(wm)
    
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
            if strcmp(wm_measure, 'ad')
                
                m_wm_LH(:, sub_count) = data_temp.ad_1(start:stop);
                sd_wm_LH(:, sub_count) = data_temp.ad_2(start:stop);
                
            elseif strcmp(wm_measure, 'fa')
                
                m_wm_LH(:, sub_count) = data_temp.fa_1(start:stop);
                sd_wm_LH(:, sub_count) = data_temp.fa_2(start:stop);
                
            elseif strcmp(wm_measure, 'md')
                
                m_wm_LH(:, sub_count) = data_temp.md_1(start:stop);
                sd_wm_LH(:, sub_count) = data_temp.md_2(start:stop);
                
            elseif strcmp(wm_measure, 'rd')
                
                m_wm_LH(:, sub_count) = data_temp.rd_1(start:stop);
                sd_wm_LH(:, sub_count) = data_temp.rd_2(start:stop);
                
            end
            
            % Grab tract name for grouping variable.
            tract_LH(:, sub_count) = {sub_contents_tractprofiles(1).name(1:end-13)};
            
            % Grab subID.
            sub_LH(sub_count) = str2num(grp_contents(i).name(end-2:end));
            
            % Grab age.
            age_LH(sub_count) = bdata.cov_age(bdata.subID == sub_LH(sub_count));
            
            clear data_temp sub_contents_tractprofiles
                            
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
            if strcmp(wm_measure, 'ad')
                
                m_wm_RH(:, sub_count) = data_temp.ad_1(start:stop);
                sd_wm_RH(:, sub_count) = data_temp.ad_2(start:stop);
                
            elseif strcmp(wm_measure, 'fa')
                
                m_wm_RH(:, sub_count) = data_temp.fa_1(start:stop);
                sd_wm_RH(:, sub_count) = data_temp.fa_2(start:stop);
                
            elseif strcmp(wm_measure, 'md')
                
                m_wm_RH(:, sub_count) = data_temp.md_1(start:stop);
                sd_wm_RH(:, sub_count) = data_temp.md_2(start:stop);
                
            elseif strcmp(wm_measure, 'rd')
                
                m_wm_RH(:, sub_count) = data_temp.rd_1(start:stop);
                sd_wm_RH(:, sub_count) = data_temp.rd_2(start:stop);
                
            end
            
            % Grab tract name for grouping variable.
            tract_RH(sub_count) = {sub_contents_tractprofiles(1).name(1:end-13)};
            
            % Grab subID.
            sub_RH(sub_count) = str2num(grp_contents(i).name(end-2:end));
            
            % Grab age group.
            age_RH(sub_count) = bdata.cov_age(bdata.subID == sub_RH(sub_count));
            
            clear data_temp sub_contents_tractprofiles
            
        end % if sub_contents empty
                                
    end % group_contents
    
    if strcmp(wm{w}, 'fa')
        
        ylim_lo = 0.10; ylim_hi = 0.70;
        
    elseif strcmp(wm{w}, 'md')
        
        ylim_lo = 0.5; ylim_hi = 1.5;
        
    end
    
    %% PLOT LEFT HEMISPHERE
    
    % Find empty cells and fill with 'empty'.
    t = find(cellfun(@isempty,tract_LH));
    tract_LH(t) = {'empty'};
    
    % Code age with color.
    max_age = max(age_LH);
    min_age = min(age_LH);
    agerange = linspace(floor(min_age), ceil(max_age), length(age_LH));
    agelength = length(agerange);
    darkblue = [0, 51, 102]/255; lightblue = [204, 229, 255]/255;
    cmap = [linspace(darkblue(1), lightblue(1), agelength)', linspace(darkblue(2), lightblue(2), agelength)', linspace(darkblue(3), lightblue(3), agelength)'];

    % Get a list of unique tract names.
    list_tract = unique(tract_LH);
    
    % Get a list of unique sub IDs.
    subID = unique(sub_LH(sub_LH~=0));
    
    if ~strcmp(list_tract{1}, 'empty')
        
        % Open a new figure and hold.
        figure(w)
        hold on;
        
        count = 0; oc_count = 0; a_count = 0;
        for s = 1:length(subID)
            
            if subID(s) ~= 0 && ~ismember(subID(s), outlier)
                
                % Find entries that are for this subject.
                s_idx = sub_LH == subID(s);
                
                % Subset the thing so that we only plot for this subject.
                t_temp = m_wm_LH(:, find(s_idx));
                
                if ~isempty(t_temp)
                    
                    count = count + 1;
                    
                    % Get index of location of this subject's age so that we
                    % can colorcode age of participant in the plot.
                    age_idx = find(agerange == interp1(agerange, agerange, age_LH(count), 'previous'));
                    
                    % Plot.
                    plot(t_temp, 'LineStyle', '-', 'LineWidth', linewidth, 'Color', [cmap(age_idx, :) alpha])
                    
                    % Collect.
                    y(:, count) = t_temp;
                    
                end
                
            end % if subID{s} ~= 0
            
        end %sub
        
    end %if strcmp
    
    % Plot means and standard deviations.
    plot(nanmean(y, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap(1, :))
    hi = nanmean(y, 2) + nanstd(y, 0, 2); lo = nanmean(y, 2) - nanstd(y, 0, 2); x = (1:size(nanmean(y, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap(1, :));
    set(hp1, 'facecolor', cmap(1, :), 'edgecolor', 'none', 'facealpha', .2);
    
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
    
    title('Left HPC to Left Parietal Cortex');
    box off;
    
    print(fullfile(rootDir, 'plots', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_LH']), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_LH']), '-depsc')
    
    hold off;
    
    clear y
    
    %% PLOT RIGHT HEMISPHERE
    
    % Find empty cells and fill with 'empty'.
    t = find(cellfun(@isempty,tract_RH));
    tract_RH(t) = {'empty'};
    
    % Code age with color.
    max_age = max(age_RH);
    min_age = min(age_RH);
    agerange = linspace(floor(min_age), ceil(max_age), length(age_LH));
    agelength = length(agerange);
    darkblue = [0, 51, 102]/255; lightblue = [204, 229, 255]/255;
    cmap = [linspace(darkblue(1), lightblue(1), agelength)', linspace(darkblue(2), lightblue(2), agelength)', linspace(darkblue(3), lightblue(3), agelength)'];
    
    % Get a list of unique tract names.
    list_tract = unique(tract_RH);
    
    % Get a list of unique sub IDs.
    subID = unique(sub_RH(sub_RH~=0));
    
    if ~strcmp(list_tract{1}, 'empty')
        
        % Open a new figure and hold.
        figure(w+2)
        hold on;
        
        count = 0;
        for s = 1:length(subID)
            
            if subID(s) ~= 0 && ~ismember(subID(s), outlier)
                
                % Find entries that are for this subject.
                s_idx = sub_RH == subID(s);
                
                % Subset the thing so that we only plot for this subject.
                t_temp = m_wm_RH(:, find(s_idx));
                
                if ~isempty(t_temp)
                    
                    count = count + 1;
                    
                    % Get index of location of this subject's age so that we
                    % can colorcode age of participant in the plot.
                    age_idx = find(agerange == interp1(agerange, agerange, age_RH(count), 'previous'));
                    
                    % Plot.
                    plot(t_temp, 'LineStyle', '-', 'LineWidth', linewidth, 'Color', [cmap(age_idx, :) alpha])
                    
                    % Collect.
                    y(:, count) = t_temp;
                    
                end
                
            end % if subID{s} ~= 0
            
        end %sub
        
    end %if strcmp
    
    % Plot means and standard deviations.
    plot(nanmean(y, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', cmap(1, :))
    hi = nanmean(y, 2) + nanstd(y, 0, 2); lo = nanmean(y, 2) - nanstd(y, 0, 2); x = (1:size(nanmean(y, 2),1))';
    hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], cmap(1, :));
    set(hp1, 'facecolor', cmap(1, :), 'edgecolor', 'none', 'facealpha', .2);
    
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
    
    title('Right HPC to Right Parietal Cortex');
    box off;
    
    print(fullfile(rootDir, 'plots', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_RH']), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{1} '_RH']), '-depsc')
    
    hold off;
    
    clear y
    
end % w












