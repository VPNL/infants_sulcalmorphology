clear; clc;
morph_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/';
addpath(morph_path)
set_participant_information

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

%% Save the processed data to CSV files

mkdir('data_csv_Figure_5b');

% Save all profiles
for r = 1:length(roi_list)
    % For each ROI, save profiles for all subjects
    roi_name = roi_list{r};
    
    % Save SD profiles
    if ~isempty(ALL_profiles_SD{r})
        SD_profiles = array2table(ALL_profiles_SD{r}, 'VariableNames', strrep(subj, '-', '_'));
        writetable(SD_profiles, ['data_csv_Figure_5b/profiles_SD_' roi_name '.csv']);
    end
    
    % Save R1 profiles
    if ~isempty(ALL_profiles_R1{r})
        R1_profiles = array2table(ALL_profiles_R1{r}, 'VariableNames', strrep(subj, '-', '_'));
        writetable(R1_profiles, ['data_csv_Figure_5b/profiles_R1_' roi_name '.csv']);
    end
    
    % Save CT profiles
    if ~isempty(ALL_profiles_CT{r})
        CT_profiles = array2table(ALL_profiles_CT{r}, 'VariableNames', strrep(subj, '-', '_'));
        writetable(CT_profiles, ['data_csv_Figure_5b/profiles_CT_' roi_name '.csv']);
    end
    
    % Save CU profiles
    if ~isempty(ALL_profiles_CU{r})
        CU_profiles = array2table(ALL_profiles_CU{r}, 'VariableNames', strrep(subj, '-', '_'));
        writetable(CU_profiles, ['data_csv_Figure_5b/profiles_CU_' roi_name '.csv']);
    end
end
