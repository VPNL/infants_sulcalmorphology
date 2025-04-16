%% Supplementary Figure 5A and C

%% 
clear; close all; clc;

save_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/Supplementary_Figure_5/';
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
    disp(['Number of ROIs to process: ', num2str(num_rois)]); 

    % Prepare arrays to store results
    Sulcus = roi_list;
    beta_Age = zeros(num_rois, 1);
    se_Age = zeros(num_rois, 1);
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
        tbl_depth_training = table(SulcalDepth(idx), rep_age(idx)', SulcalSpan(idx),...
           Thickness(idx), Curvature(idx),R1_gray(idx),rep_group(idx)',...
           'VariableNames', {'SD', 'Age', 'SP', 'CT', 'CU', 'R1', 'Baby'});


        % Fit the linear mixed-effects model
        model = fitlme(tbl_depth_training, 'SD ~ 1 + (Age + SP + CT + CU + R1) + (1|Baby)');
        disp(['Model fit for ROI: ', currentROI]); 

        % Store beta and standard error values
        beta_Age(roi_idx) = model.Coefficients.Estimate(2); 
        se_Age(roi_idx) = model.Coefficients.SE(2);         
        beta_SP(roi_idx) = model.Coefficients.Estimate(3);  
        se_SP(roi_idx) = model.Coefficients.SE(3);          
        beta_CT(roi_idx) = model.Coefficients.Estimate(4);  
        se_CT(roi_idx) = model.Coefficients.SE(4);          
        beta_CU(roi_idx) = model.Coefficients.Estimate(5);  
        se_CU(roi_idx) = model.Coefficients.SE(5);          
        beta_R1(roi_idx) = model.Coefficients.Estimate(6);  
        se_R1(roi_idx) = model.Coefficients.SE(6);          
    end

    % Create the table with the specified columns
    results_table = table(Sulcus', beta_Age, se_Age, beta_SP, se_SP, beta_CT, se_CT, ...
                      beta_CU, se_CU, beta_R1, se_R1, ...
                      'VariableNames', {'Sulcus', 'Beta_Age', 'SE_Age', 'Beta_SP', 'SE_SP', ...
                                        'Beta_CT', 'SE_CT', 'Beta_CU', 'SE_CU','Beta_R1', 'SE_R1'});

    % Save the table to the specified path
    output_file = fullfile(save_path, [hemisphere '_Figure_4_Normalized_Betas_age_added.xlsx']);
    writetable(results_table, output_file);



    % Plot the beta values
    ci = 1.96; % For 95% confidence intervals
    fig = figure('Color', 'white'); 
    set(fig, 'Position', [343, 168, 3328, 703]);
    loc = 0.05;

    % Prepare beta and SE values for plotting
    predictors = {'Age', 'SP', 'CT', 'CU','R1'};
    beta_values = [beta_Age, beta_SP, beta_CT, beta_CU, beta_R1];
    se_values   = [se_Age,  se_SP,  se_CT,  se_CU,  se_R1];


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
        ylim([-0.7 1])

        % Customize appearance
        if roi_idx == 1
            set(gcf, 'Color', 'w');
            box off;
            set(gca, 'XTick', []);
            yticklabels("");
        else
            set(gcf, 'Color', 'w');
            box off;
            %set(gca, 'XTick', []);
            set(gca, 'YColor', 'none'); 
        end

        hold off;
        loc = loc + 0.05;
    end

    % % Save the figure
    fig_output_file = fullfile(save_path, [hemisphere '_ROI_Normalized_BetaPlot_age_added.png']);
    saveas(fig, fig_output_file);
end
