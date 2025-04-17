%% FIGURE 4C

% Note: All parameters have been normalized prior to model fitting and visualization 
% to allow for direct comparison across metrics, which are originally measured 
% in different units and ranges. 

%% 
clear; close all; clc;

save_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/Figure_4c/';
set_participant_information

hemispheres = {'LH', 'RH'}; 
rep_age = repmat(logage,1,12);
rep_group = repmat(group,1,12);

for h = 1:length(hemispheres)
    hemisphere = hemispheres{h};

    % Load hemisphere-specific data
    data_file = [hemisphere '_all_parameters_normalized.csv'];
    dataTbl = readtable(data_file);
    SulcalDepth = dataTbl.SulcalDepth;
    SulcalSpan = dataTbl.SulcalWidth;
    Curvature = dataTbl.Curvature;
    R1_gray = dataTbl.R1_gray;
    Thickness = dataTbl.Thickness;
    ROI = dataTbl.ROI;

    num_rois = length(roi_list);

    % Prepare arrays to store results
    Sulcus = roi_list;
    beta_SP = zeros(num_rois, 1);
    se_SP = zeros(num_rois, 1);
    beta_CT = zeros(num_rois, 1);
    se_CT = zeros(num_rois, 1);
    beta_CU = zeros(num_rois, 1);
    se_CU = zeros(num_rois, 1);
    beta_R1 = zeros(num_rois, 1);
    se_R1 = zeros(num_rois, 1);
    

    % Loop through each ROI and fit the model, storing the beta and SE values
    for roi_idx = 1:num_rois
        currentROI = roi_list{roi_idx};
        idx = strcmp(ROI, currentROI);

        % Create the table for model fitting
        tbl_depth_training = table(SulcalDepth(idx), SulcalSpan(idx), Thickness(idx), ...
                                   Curvature(idx), R1_gray(idx), rep_age(idx)', rep_group(idx)',...
                                   'VariableNames', {'SD', 'SP', 'CT', 'CU', 'R1', 'Age', 'Baby'});

        % Fit the linear mixed-effects model
        model = fitlme(tbl_depth_training, 'SD ~ 1 + (SP + CT + CU + R1) + (1|Baby)');
        disp(['Model fit for ROI: ', currentROI]); 

        % Store beta and standard error values
        beta_SP(roi_idx) = model.Coefficients.Estimate(2); % SP beta
        se_SP(roi_idx) = model.Coefficients.SE(2);         % SP SE
        beta_CT(roi_idx) = model.Coefficients.Estimate(2); % CT beta
        se_CT(roi_idx) = model.Coefficients.SE(2);         % CT SE
        beta_CU(roi_idx) = model.Coefficients.Estimate(3); % CU beta
        se_CU(roi_idx) = model.Coefficients.SE(3);         % CU SE
        beta_R1(roi_idx) = model.Coefficients.Estimate(4); % R1 beta
        se_R1(roi_idx) = model.Coefficients.SE(4);         % R1 SE
    end

    % Create the table with the specified columns
    results_table = table(Sulcus', beta_CT, se_CT,  beta_CU, se_CU, beta_R1, se_R1,...
                          'VariableNames', {'Sulcus', 'Beta_CT', 'SE_CT', 'Beta_CU', 'SE_CU','Beta_R1', 'SE_R1'});

    % Save the table to the specified path
    output_file = fullfile(save_path, [hemisphere '_Figure_4_Normalized_Betas.xlsx']);
    writetable(results_table, output_file);

    % Plot the beta values
    ci = 1.96; % For 95% confidence intervals
    predictors = {'CT', 'CU','R1'};
    fig = figure('Color', 'white'); 
    set(fig, 'Position', [343, 168, 3328, 703]);
    loc = 0.05;

    % Prepare beta and SE values for plotting
    beta_values = [beta_CT, beta_CU, beta_R1];
    se_values = [se_CT, se_CU, se_R1];

    for roi_idx = 1:num_rois
        aa = subplot(1, num_rois, roi_idx); hold;
        pos = get(aa, 'Position');
        pos(1) = loc;
        set(aa, 'Position', pos);

        % Plot bar chart
        bar(beta_values(roi_idx, :), 'FaceColor', color(roi_idx, :), 'EdgeColor', color(roi_idx, :));
        hold on;
        % Add error bars
        errorbar(1:length(predictors), beta_values(roi_idx, :), ...
                 se_values(roi_idx, :), ... 
                 'Color', [0.6 0.6 0.6], 'LineStyle', 'none', 'LineWidth', 1);

        set(gca, 'XTick', [])
        ylim([-0.6 1])

        % Customize appearance
        if roi_idx == 1
            set(gcf, 'Color', 'w');
            box off;
            set(gca, 'XTick', []);
            yticklabels("");
        else
            set(gcf, 'Color', 'w');
            box off;
            set(gca, 'XTick', []);
            set(gca, 'YColor', 'none'); 
        end

        hold off;
        loc = loc + 0.05;
    end

    %Save the figure
    fig_output_file = fullfile(save_path, [hemisphere '_ROI_Normalized_BetaPlot.png']);
    saveas(fig, fig_output_file);
end
