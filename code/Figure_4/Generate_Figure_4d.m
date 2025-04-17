%% Figure 4D

%%
clear; clc;

morph_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/';
addpath(morph_path)
set_participant_information

roi_list = {'calcarine','pos','insula','central','cos','sts','sfs','ips','loc','ifs','ots_lat','its'};
save_path = fullfile(morph_path, 'Figure_4d');

hemispheres = {'lh', 'rh'};

for h = 1:length(hemispheres)
    hemi = hemispheres{h};

    % === Load hemisphere-specific data from CSV ===
    SD_tbl = readtable([morph_path 'SD_' hemi '.csv']);
    R1gray_tbl = readtable([morph_path 'R1gray_' hemi '.csv']);
    thickness_tbl = readtable([morph_path 'thickness_' hemi '.csv']);
    SA_tbl = readtable([morph_path 'SP_' hemi '.csv']);
    curvature_tbl = readtable([morph_path 'curvature_' hemi '.csv']);

    % === Flatten tables into vectors (subjects x ROIs → 1 x (subjects*ROIs)) ===
    SD = reshape(table2array(SD_tbl), 1, []);
    R1gray = reshape(table2array(R1gray_tbl), 1, []);
    thickness_data = reshape(table2array(thickness_tbl), 1, []);
    SP = reshape(table2array(SA_tbl), 1, []);
    curvature = reshape(table2array(curvature_tbl), 1, []);

    % === Repeat labels for each ROI ===
    ROI_id = 1:length(roi_list);
    repeated_ROI = repelem(ROI_id, 79);

    ROI_emergence = [ones(1,79)*16 ones(1,79)*16 ones(1,79)*18 ones(1,79)*20 ones(1,79)*20 ...
                     ones(1,79)*24 ones(1,79)*24 ones(1,79)*24 ones(1,79)*24 ...
                     ones(1,79)*28 ones(1,79)*28 ones(1,79)*28];

    rep_age = repmat(logage, 1, length(roi_list));
    rep_group = repmat(group, 1, length(roi_list));
    frame1 = repmat(subj, 1, length(roi_list));

    % === Build table ===
    tbl_depth_training = table(rep_age', R1gray', thickness_data', SD', SP', curvature', ...
        repeated_ROI', ROI_emergence', rep_group', frame1', ...
        'VariableNames', {'Age','R1','CT','SD','SP','CU','ROI','Emergence','Baby','Bname'});

    % === Initialize output containers ===
    RMSE_results = zeros(length(roi_list), 1);
    STD_results = zeros(length(roi_list), 1);
    beta_values = cell(length(roi_list), 1);
    final_predictions = cell(length(roi_list), 1);
    actual_depth = cell(length(roi_list), 1);

    % === Model loop ===
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
        RMSE_results(roi_idx) = sqrt(mean(errors.^2));
        STD_results(roi_idx) = std(errors);
    end

    % === Plot predicted vs actual ===
    fig = figure('Color','white');
    set(fig, 'position', [343 168 3328 703]);
    loc = 0.1;

    for roi_idx = 1:length(roi_list)
        aa = subplot(1,12, roi_idx); hold;
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

    % === Save figure ===
    full_file_path = fullfile(save_path, [hemi '_roi_predicted_observed_with_age.tiff']);
    saveas(gcf, full_file_path);

    % === Print results ===
    fprintf('\n---- Results for Hemisphere: %s ----\n', hemi);
    for roi_idx = 1:length(roi_list)
        fprintf('ROI: %s - RMSE: %.4f ± %.4f\n', roi_list{roi_idx}, RMSE_results(roi_idx), STD_results(roi_idx));
        fprintf('Beta values for %s: ', roi_list{roi_idx});
        disp(beta_values{roi_idx}');
    end
end
