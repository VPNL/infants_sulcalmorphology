clear; clc;

set_participant_information

% Set hemisphere to analyze
hemi = 'lh'; % Change to 'rh' for right hemisphere


% Load the pre-processed beta values
betable = readtable(['/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/data_csv_Figure_5c/', hemi, '_normalized_betas.csv']);
    
% Create Figure 5c
fig = figure('Color', 'white'); 
set(fig, 'Position', [343, 168, 3328, 703]);
loc = 0.05;
count = 1;
num_rois = length(roi_list);

for roi_idx = 1:num_rois
    aa = subplot(1, num_rois, roi_idx); hold;
    pos = get(aa, 'Position');
    pos(1) = loc;
    set(aa, 'Position', pos);
    
    start = (roi_idx - 1) * 3 + 1;
    ende = roi_idx * 3;

    b = bar(betable.EST(start:ende));
    hold on;
    errorbar(1:3, betable.EST(start:ende), betable.SE(start:ende), betable.SE(start:ende), 'linestyle', 'none', 'color', [.6 .6 .6]);
    b.FaceColor = color(count,:); b.LineWidth = 1.5;
    b.EdgeColor = color(count,:);
    count = count + 1;
    set(gcf, 'Color', 'w');
    ylim([-0.8 1.3])

    % Customize appearance to make figure-ready 
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
    normAxes
    loc = loc + 0.05; 
end