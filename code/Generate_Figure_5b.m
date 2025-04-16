%% Figure 5b

% The script has two parts:
% 1. REFERENCE CODE (not runnable without own data): Shows how the original profiles 
%    were created from NIFTI files. This section requires access to the original 
%    brain scan data and is included for reference only.

% 2. RUNNABLE CODE: Creates Figure 5b using pre-processed data saved in 
%    'data_csv_Figure_5b' directory. 

%% REFERENCE CODE
clear; clc;
morph_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/';
addpath(morph_path)
set_participants_79

hemi= 'lh';

Correl_R1=[]; Correl_CT=[]; Correl_CU=[];

R1_data=[]; CU_data=[]; SD_data=[]; CT_data=[];
ALL_profiles_SD = []; ALL_profiles_CT = []; ALL_profiles_CU = []; ALL_profiles_R1 = [];

% Load data for each subject 
for i= 1:length(subj) 
    
    % Get subject info 
    whichmonth_path(i) = extractBetween(subj{i}, '_', '_');
    whichbaby{i} = extractBefore(subj{i}, '_');

    % Load R1 data
    cd(['/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/' whichbaby{i} '/' whichmonth_path{i} '/preprocessed_aligned/qmri/'])
    file = ([hemi,'.SEIR_R1_aligned_1mm_fsaverage.nii.gz']);
    tempfile = readFileNifti(file);
    data = reshape((tempfile.data) ,1 , 256*256*256);
    R1_data(:,i) = data;
    
    % Load sulcal depth data 
    cd(['/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains/' subj{i} '/surf'])
    file = [subj{i} '_sulc_fsaverage_',hemi,'.nii.gz'];
    tempfile = readFileNifti(file);
    data = reshape(tempfile.data ,1 , 256*256*256);
    SD_data(:,i) = data; 
    
    % Load thickness data
    file = ([subj{i} '_THICK_fsaverage_',hemi,'.nii.gz']);
    tempfile = readFileNifti(file);
    data = reshape(tempfile.data ,1 , 256*256*256);
    CT_data(:,i) = data;
    
    % Load curvature data
    file = ([subj{i} '_CURV_fsaverage_',hemi,'.nii.gz']);
    tempfile = readFileNifti(file);
    data = reshape(tempfile.data ,1 , 256*256*256);
    CU_data(:,i) = data;    
end

% Create profiles for each ROI 
for r =  1:length(roi_list)
    % get ROI labels 
    roi_name = roi_list{r};
    ROIvals = readFileNifti([FSdir '/fsAverage/label/sarah_coarse_labels/', hemi,'.', roi_name '.coarse.st.label.nii.gz']);
    ROIlabel = reshape(ROIvals.data ,1 , 256*256*256);
    
    SD_tempdata =[];  %% this is the SD profile
    R1_tempdata =[];  %% this is the R1 profile
    CT_tempdata = []; %% this is for CT profile
    CU_tempdata= [];  %% this is for CU profile
    
    % Extract data for this ROI 
    CU_tempdata = CU_data(find(ROIlabel==1),:);
    CT_tempdata = CT_data(find(ROIlabel==1),:);
    SD_tempdata = SD_data(find(ROIlabel==1),:);
    R1_tempdata = R1_data(find(ROIlabel==1),:);
    
    % Clean up extreme values 
    for s=1:length(subj)
        temp= R1_tempdata(:,s);
        outlier_vals= find(temp>2);
        R1_tempdata(outlier_vals,s)= NaN;

        temp= CU_tempdata(:,s);
        outlier_vals= find(temp<-.1);
        CU_tempdata(outlier_vals,s)= NaN;      
    end
    
    % Get 3D dimenstions
    ROIlabel = reshape(ROIlabel ,256,256,256);
    [ii jj kk]=ind2sub(size(ROIlabel), find(ROIlabel==1)); %
    
    % Create profiles based on ROI direction
    if r==2 % Superior to inferior
        profile_SD = []; profile_R1 = []; profile_CT = []; profile_CU = [];
        
        for ll = min(jj):1:max(jj)
            profile_SD = [profile_SD; nanmean(SD_tempdata(find(jj==ll),:),1)];
            profile_R1 = [profile_R1; nanmean(R1_tempdata(find(jj==ll),:),1)];
            profile_CT = [profile_CT; nanmean(CT_tempdata(find(jj==ll),:),1)];
            profile_CU = [profile_CU; nanmean(CU_tempdata(find(jj==ll),:),1)];
        end
    
        
    elseif  r == 4 || r == 8  %% Left to right 

        profile_SD = []; profile_R1 = []; profile_CT = []; profile_CU = [];
        for ll = min(ii):1:max(ii)
            profile_SD = [profile_SD; nanmean(SD_tempdata(find(ii==ll),:),1)];
            profile_R1 = [profile_R1; nanmean(R1_tempdata(find(ii==ll),:),1)];
            profile_CT = [profile_CT; nanmean(CT_tempdata(find(ii==ll),:),1)];
            profile_CU = [profile_CU; nanmean(CU_tempdata(find(ii==ll),:),1)];
        end
    else % Anterior to posterior
        profile_SD = []; profile_R1 = []; profile_CT = []; profile_CU = [];
        for ll = min(kk):1:max(kk)
            profile_SD = [profile_SD; nanmean(SD_tempdata(find(kk==ll),:),1)];
            profile_R1 = [profile_R1; nanmean(R1_tempdata(find(kk==ll),:),1)];
            profile_CT = [profile_CT; nanmean(CT_tempdata(find(kk==ll),:),1)];
            profile_CU = [profile_CU; nanmean(CU_tempdata(find(kk==ll),:),1)];
        end
    end

    % Store profiles
    ALL_profiles_SD{r} = profile_SD;
    ALL_profiles_R1{r} = profile_R1;
    ALL_profiles_CU{r} = profile_CU;
    ALL_profiles_CT{r} = profile_CT;

end 

%% RUNNABLE CODE 

r = 6; % STS ROI
i = 55; % Specific subject index for figure

% ROI name
roi_name = roi_list{r};

% Load profiles for this ROI
SD_profiles = readtable(['data_csv_Figure_5b/profiles_SD_' roi_name '.csv']);
R1_profiles = readtable(['data_csv_Figure_5b/profiles_R1_' roi_name '.csv']);
CT_profiles = readtable(['data_csv_Figure_5b/profiles_CT_' roi_name '.csv']);
CU_profiles = readtable(['data_csv_Figure_5b/profiles_CU_' roi_name '.csv']);

% Convert tables to arrays
SD_profile = table2array(SD_profiles(:, i));
R1_profile = table2array(R1_profiles(:, i));
CT_profile = table2array(CT_profiles(:, i));
CU_profile = table2array(CU_profiles(:, i));

% Create Figure 5b
figure;
% SD vs CU
subplot(3,1,1); hold;
plot(normalize(SD_profile), 'k', 'linewidth', 2);
plot(normalize(CU_profile), ':', 'Color', [0 128 255]/255, 'linewidth', 2);
[r1, p] = corrcoef(SD_profile, CU_profile, 'rows', 'complete');
title([char(subj(i)), ' R: ', num2str(r1(1,2)), ' p: ', num2str(p(1,2))]);
hold off;

% SD vs R1
subplot(3,1,2); hold;
plot(normalize(SD_profile), 'k', 'linewidth', 2);
plot(normalize(R1_profile), ':', 'Color', [102 204 0]/255, 'linewidth', 2);
[r1, p] = corrcoef(SD_profile, R1_profile, 'rows', 'complete');
title(['R: ', num2str(r1(1,2)), ' p: ', num2str(p(1,2))]);
hold off;

% SD vs CT
subplot(3,1,3); hold;
plot(normalize(SD_profile), 'k', 'linewidth', 2);
plot(normalize(CT_profile), ':', 'Color', [218 84 0]/255, 'linewidth', 2);
[r1, p] = corrcoef(SD_profile, CT_profile, 'rows', 'complete');
title(['R: ', num2str(r1(1,2)), ' p: ', num2str(p(1,2))]);
hold off;