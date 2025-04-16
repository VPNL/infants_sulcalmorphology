%this script compares the overlap between ROIs drawn by hand and freesurfer
%projected ROIs from the adult average brain 

clear all; clc;

% 
subj = {
    'bb02_mri0_mask', 'bb11_mri0_mask', 'bb56_mri0_mask', 'bb22_mri0_mask',...
    'bb62_mri0_mask', 'bb07_mri0_mask', 'bb67_mri0_mask', 'bb28_mri0_mask',...
    'bb14_mri0_mask', 'bb37_mri0_mask', 'bb15_mri3_mask', 'bb11_mri3_mask',...
    'bb02_mri3_mask', 'bb40_mri3_mask', 'bb75_mri3_mask', 'bb24_mri3_mask',...
    'bb28_mri3_mask', 'bb07_mri3_mask', 'bb30_mri3_mask', 'bb52_mri3_mask',...
    'bb40_mri6_mask', 'bb66_mri6_mask', 'bb32_mri6_mask', 'bb82_mri6_mask',...
    'bb19_mri6_mask', 'bb47_mri6_mask', 'bb11_mri6_mask', 'bb12_mri6_mask',...
    'bb62_mri6_mask', 'bb05_mri6_mask', 'bb28_mri12_mask', 'bb56_mri12_mask',...
    'bb62_mri12_mask', 'bb55_mri12_mask', 'bb37_mri12_mask', 'bb32_mri12_mask',...
    'bb40_mri12_mask', 'bb61_mri12_mask', 'bb42_mri12_mask', 'bb47_mri12_mask'
};


age = [29, 24, 9, 30, 44, 37, 21, 26, 31, 50, 104, 78, 85, 123, 106, 94, 94, 95, 82, 105, 190, 201, 182, 198, 177, 170, 167, 181, 182, 189, 368, 418, 392, 398, 373, 360, 372, 385, 399, 393];

color = [255 102 102; 255 178 102; 255 255 102; 153 255 51; 51 255 51; 51 255 153; 51 255 255; 51 153 255; 51 51 255; 153 51 255; 255 51 255; 255 51 153]/255;

roi_list = {'calcarine','pos','insula','central','cos','sts','sfs','ips','loc','ifs','ots','its'};

hemisphere = {'lh','rh'};

ROI_type = {'coarse','groundtruth'};

i = 0;

k = 0;


for s = 1:length(subj)

    src_folder = ['/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains/' subj{s} '/label/sarah_coarse/'];
    dest_folder = ['/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains/' subj{s} '/label/'];

    for h = 1:length(hemisphere)
        for r = 12%1:length(roi_list)
            for t = 1:length(ROI_type)
                if strcmpi(ROI_type{t}, 'groundtruth')
                    fname = [hemisphere{h} '.' roi_list{r} '.' ROI_type{t} '.label'];
                else
                    fname = [hemisphere{h} '.' roi_list{r} '.' ROI_type{t} '.st.label'];
                end

                % Construct full paths for source and destination
                src_file = fullfile(src_folder, fname);
                dest_file = fullfile(dest_folder, fname);
                
                % Check if file exists before trying to move
                if isfile(src_file)
                    movefile(src_file, dest_file);
                else
                    fprintf('File %s does not exist for subject %s.\n', fname, subj{s});
                end
            end
        end
    end
end 

%% 


      
for t = 1:length(ROI_type)
        for s = 1:length(subj)
            cd(['/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains/' subj{s} '/label/sarah_coarse/'])
            for h = 1:length(hemisphere)
                for r = 1:length(roi_list)
                    if strcmpi(ROI_type{t}, 'groundtruth')
                        lname = ([hemisphere{h} '.' roi_list{r} '.' ROI_type{t}]);
                    else
                        lname = ([hemisphere{h} '.' roi_list{r} '.' ROI_type{t} '.st']);
                    end
                    sname = (subj{s});
                    data = read_label(sname, lname);
                    if strcmpi(ROI_type{t}, 'groundtruth')
                        i = i + 1;
                        data_groundtruth{i} = data(:, 1);
                    elseif strcmpi(ROI_type{t}, 'coarse')
                        k = k + 1;
                        data_coarse{k} = data(:, 1);
                    end 
                end
            end
        end
    end


 
    %% 

for l = 1:length(data_coarse)
    disp(['Size of data_groundtruth: ', num2str(length(data_groundtruth))]);
disp(['Size of data_coarse: ', num2str(length(data_coarse))]);
disp(['Current iteration (l): ', num2str(l)]);

    overlap_coarse = intersect(data_groundtruth{l}, data_coarse{l});
    DC_coarse(l) = 2 * (length(overlap_coarse)) / (length(data_groundtruth{l}) + length(data_coarse{l}));
end 

for r = 1:length(roi_list)
    D_coarse_lh(r,:) = DC_coarse(r:length(roi_list)*length(hemisphere):end);
    DC_coarse_mean_lh(r) = mean(D_coarse_lh(r,:));
  
    D_coarse_rh(r,:) = DC_coarse(r+length(roi_list):length(roi_list)*length(hemisphere):end);
    DC_coarse_mean_rh(r) = mean(D_coarse_rh(r,:));
  
    DC_coarse_sem_lh(r) = std(D_coarse_lh(r,:)) / sqrt(length(D_coarse_lh(r,:)));
    DC_coarse_sem_rh(r) = std(D_coarse_rh(r,:)) / sqrt(length(D_coarse_rh(r,:)));
end

figure;
set(gcf, 'color', 'white');

subplot(2,1,1);
x = 1:length(roi_list);

y = [DC_coarse_mean_rh' DC_coarse_mean_lh'];

lighter_color = color + (1-color) * 0.5;  % Mix with white to get a lighter shade

% Stack colors for rh and lh in a single matrix
all_colors = [color; lighter_color];

hbar = bar(x, y, 'grouped');

for k1 = 1:size(y,2) 
    for k2 = 1:size(y,1) 
        hbar(k1).FaceColor = 'flat';  
        if k1 == 1 % Left hemisphere
            hbar(k1).CData(k2, :) = color(k2, :); 
        else % Right hemisphere
            hbar(k1).CData(k2, :) = lighter_color(k2, :);
        end
    end
end


err = [DC_coarse_sem_rh' DC_coarse_sem_lh'];

hold on;

barWidth = 0.4; 
ngroups = size(y, 1);
nbars = size(y, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    % Position of bar centers, relative to each error bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:, i), err(:, i), '.k');
end

box off;
set(gca, 'XTickLabel', []);  % Removes the x-labels
set(gca, 'XTick', []);      % Removes the x-tick marks
set(gca, 'YTickLabel', []);  % Removes the x-labels
set(gca, 'YTick', []);  

%% 
subjects_per_age_group = 10;
age_groups = [0, 3, 6, 12];
num_sulci = 12; % Number of ROIs/sulci
num_subjects = 40;
hemispheres = {'Left', 'Right'};
num_hemispheres = numel(hemispheres);

% Preallocate arrays for ANOVA data
age_group_data = zeros(num_sulci * num_subjects * num_hemispheres, 1);
hemisphere_data = cell(num_sulci * num_subjects * num_hemispheres, 1);
sulcus_data = cell(num_sulci * num_subjects * num_hemispheres, 1);
dice_coefficients = zeros(num_sulci * num_subjects * num_hemispheres, 1);

% Fill in the data
idx = 1;
for s = 1:num_sulci
    for subj = 1:num_subjects
        for h = 1:num_hemispheres
            age_group_data(idx) = age_groups(ceil(subj / subjects_per_age_group));
            hemisphere_data{idx} = hemispheres{h};
            sulcus_data{idx} = sprintf('Sulcus_%d', s);
            if h == 1 % Left hemisphere
                dice_coefficients(idx) = D_coarse_lh(s, subj);
            else % Right hemisphere
                dice_coefficients(idx) = D_coarse_rh(s, subj);
            end
            idx = idx + 1;
        end
    end
end
% Convert numerical age group data to categorical for ANOVA
age_group_data_categorical = categorical(age_group_data);

% Perform 3-way ANOVA
[p, tbl, stats, terms] = anovan(dice_coefficients, {age_group_data_categorical, hemisphere_data, sulcus_data},...
    'model', 'interaction', 'varnames', {'Age', 'Hemisphere', 'Sulcus'});

%% 

% Initialize a structure to store results
roi_stats = struct();

% Loop over each ROI and compute statistics
for r = 1:length(roi_list)
    % Left hemisphere
    mean_lh = mean(D_coarse_lh(r, :));
    std_lh = std(D_coarse_lh(r, :));
    roi_stats.lh.(roi_list{r}) = sprintf('%.2f: %.2f', mean_lh, std_lh);

    % Right hemisphere
    mean_rh = mean(D_coarse_rh(r, :));
    std_rh = std(D_coarse_rh(r, :));
    roi_stats.rh.(roi_list{r}) = sprintf('%.2f: %.2f', mean_rh, std_rh);
end

% Compute mean and standard deviation for the left hemisphere
mean_dice_lh = mean(D_coarse_lh(:), 'omitnan'); % Mean of LH Dice coefficients
std_dice_lh = std(D_coarse_lh(:), 'omitnan');   % Standard deviation of LH Dice coefficients

% Compute mean and standard deviation for the right hemisphere
mean_dice_rh = mean(D_coarse_rh(:), 'omitnan'); % Mean of RH Dice coefficients
std_dice_rh = std(D_coarse_rh(:), 'omitnan');   % Standard deviation of RH Dice coefficients

% Display the results
disp(['Mean Dice Coefficient for Left Hemisphere: ', num2str(mean_dice_lh, '%.4f')]);
disp(['Standard Deviation of Dice Coefficients for Left Hemisphere: ', num2str(std_dice_lh, '%.4f')]);
disp(['Mean Dice Coefficient for Right Hemisphere: ', num2str(mean_dice_rh, '%.4f')]);
disp(['Standard Deviation of Dice Coefficients for Right Hemisphere: ', num2str(std_dice_rh, '%.4f')]);
