% Supplementary Figure 2B
clear all; clc;

data_dir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/roi_labels_dice';

set_participant_information

% Load mean and SEM data
dice_means = readtable(fullfile(data_dir, 'dice_means_sem.csv'));

% Extract means and SEMs
DC_coarse_mean_lh = dice_means.Mean_LH;
DC_coarse_mean_rh = dice_means.Mean_RH;
DC_coarse_sem_lh = dice_means.SEM_LH;
DC_coarse_sem_rh = dice_means.SEM_RH;

% Create the figure
figure;
set(gcf, 'color', 'white');
set(gcf, 'Position', [100, 100, 1200, 400]); 

% Create the bar plot
x = 1:length(roi_list);
y = [DC_coarse_mean_rh DC_coarse_mean_lh];
err = [DC_coarse_sem_rh DC_coarse_sem_lh];

% Calculate lighter color for right hemisphere
lighter_color = color + (1-color) * 0.5;  % Mix with white to get a lighter shade

% Create the grouped bar chart
hbar = bar(x, y, 'grouped');

% Set bar colors
for k1 = 1:size(y,2) 
    for k2 = 1:size(y,1) 
        hbar(k1).FaceColor = 'flat';  
        if k1 == 1 % Right hemisphere
            hbar(k1).CData(k2, :) = color(k2, :); 
        else % Left hemisphere
            hbar(k1).CData(k2, :) = lighter_color(k2, :);
        end
    end
end

% Add error bars
hold on;
barWidth = 0.4; 
ngroups = size(y, 1);
nbars = size(y, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    x_pos = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x_pos, y(:, i), err(:, i), '.k');
end

% Format
set(gca, 'XTick', 1:length(roi_list), 'XTickLabel', roi_list, 'FontSize', 10);
set(gca, 'XTickLabelRotation', 45); % Rotate labels for better readability

legend([hbar(1), hbar(2)], {'Right Hemisphere', 'Left Hemisphere'}, 'Location', 'best');

grid off;
box off;

% Adjust y-axis to better show the data
ylim([0 max(max(y + err)) * 1.1]);
