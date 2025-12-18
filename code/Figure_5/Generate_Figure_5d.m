%% Generate Figure 5d
% Note: This model was trained and evaluated using raw (non-normalized) values 
% to preserve interpretability in millimeters

%%
clear; clc;

set_participant_information
addpath(morph_path)

% Paths to change
save_path = fullfile(morph_path, 'REVISION','Figure_4d');
file_path = fullfile(morph_path, 'REVISION');

hemispheres = {'lh', 'rh'};

for h = 1:length(hemispheres)
    hemi = hemispheres{h};

    % Load data 
    SD_tbl = readtable([file_path filesep 'SULC_' hemi '.csv']);
    R1_tbl = readtable([file_path filesep 'R1_' hemi '.csv']);
    CT_tbl = readtable([file_path filesep 'THICK_' hemi '.csv']);
    SP_tbl = readtable([file_path filesep 'AREA_' hemi '.csv']);
    CU_tbl = readtable([file_path filesep 'CURV_' hemi '.csv']);

    % Convert to array and skip first column (first column is just the
    % subject name)
    SD = reshape(table2array(SD_tbl(:,2:end)), 1, []);  
    R1 = reshape(table2array(R1_tbl(:,2:end)), 1, []);  
    CT = reshape(table2array(CT_tbl(:,2:end)), 1, []);  
    SP = reshape(table2array(SP_tbl(:,2:end)), 1, []);  
    CU = reshape(table2array(CU_tbl(:,2:end)), 1, []);  

    % Repeat labels for each ROI
    ROI_id = 1:length(roi_list);
    repeated_ROI = repelem(ROI_id, 79);

    rep_age = repmat(logage, 1, length(roi_list));
    rep_group = repmat(group, 1, length(roi_list));
    frame1 = repmat(subj, 1, length(roi_list));

    % Build table for training
    tbl_depth_training = table(rep_age', R1', CT', SD', SP', CU', ...
        repeated_ROI', rep_group', frame1', ...
        'VariableNames', {'Age','R1','CT','SD','SP','CU','ROI','Baby','Bname'});

    % Initialize results
    RMSE_results = zeros(length(roi_list), 1);
    STD_results = zeros(length(roi_list), 1);
    beta_values = cell(length(roi_list), 1);
    final_predictions = cell(length(roi_list), 1);
    actual_depth = cell(length(roi_list), 1);

    for roi_idx = 1:length(roi_list)
        roi_name = roi_list{roi_idx};
        roi_color = color(roi_idx, :);

        roi_data = tbl_depth_training(tbl_depth_training.ROI == roi_idx, :);
        predictions = zeros(size(roi_data, 1), 1);
        model_coefficients = [];

        for i = 1:size(roi_data, 1)
            test_data = roi_data(i, :);
            train_data = roi_data(~strcmp(roi_data.Bname, test_data.Bname{1}), :);

            model = fitlme(train_data, 'SD ~ 1 + (SP + R1 + CU + CT) + (1|Baby)');
            predictions(i) = predict(model, test_data);

            if isempty(model_coefficients)
                model_coefficients = fixedEffects(model);
            end
        end

        beta_values{roi_idx} = model_coefficients;
        final_predictions{roi_idx} = predictions;
        actual_depth{roi_idx} = roi_data.SD;

        errors = predictions - roi_data.SD;
        mean_error = mean(errors)
        fprintf('Mean error for %s: %.4f\n', roi_list{roi_idx}, mean_error);
        histogram(errors);
        title(sprintf('Error Distribution for %s', roi_list{roi_idx}));
        xlabel('Error');
        ylabel('Frequency');
        direct_RMSE = sqrt(mean(errors.^2));
        fprintf('Direct RMSE for %s: %.4f\n', roi_list{roi_idx}, direct_RMSE);
    
        mean_prediction = mean(predictions);
        mean_actual = mean(roi_data.SD);
        fprintf('Mean prediction: %.4f, Mean actual: %.4f\n', mean_prediction, mean_actual);
    
        figure;
        scatter(roi_data.SD, predictions);
        hold on;
        plot([min(roi_data.SD), max(roi_data.SD)], [min(roi_data.SD), max(roi_data.SD)], 'r--');
        xlabel('Actual');
        ylabel('Predicted');
        title(sprintf('Predictions vs Actuals for %s', roi_list{roi_idx}));
  
        RMSE_results(roi_idx) = sqrt(mean(errors.^2));
        STD_results(roi_idx) = std(errors);
    end
    
    % Plot predicted vs actual 
    fig = figure('Color','white');
    set(fig, 'position', [343 168 3328 703]);
    loc = 0.1;

    for roi_idx = 1:length(roi_list)
        aa = subplot(1,length(roi_list), roi_idx); hold;
        pos = get(aa,'Position'); pos(1)=loc; set(aa,'Position',pos);

        scatter(actual_depth{roi_idx}, final_predictions{roi_idx}, ...
            70, 'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerFaceColor', color(roi_idx,:), ...
            'DisplayName', roi_list{roi_idx});

        x = min(actual_depth{roi_idx}):0.1:max(actual_depth{roi_idx});
        plot(x, x, 'k','LineWidth', 2); % identity line

        axis([0 8 0 8]); axis square; box off;
        if roi_idx ~= 1
            set(gca,'YColor','none');
        end
        loc = loc + 0.055;
        hold off;
    end

    % Save figure 
    full_file_path = fullfile(save_path, [hemi '_roi_predicted_observed_with_age.png']);
    saveas(gcf, full_file_path);
end