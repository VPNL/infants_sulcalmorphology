%% Figure 5c

clear
morph_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/';
addpath(morph_path)
cd(morph_path)
set_participant_information

hemi = 'rh';
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
for r = 1:length(roi_list)
    % Get ROI labels indices
    roi_name = roi_list{r};
    ROIvals = readFileNifti([FSdir '/fsAverage/label/sarah_coarse_labels/', hemi,'.', roi_name '.coarse.st.label.nii.gz']);
    ROIlabel = reshape(ROIvals.data ,1 , 256*256*256);
    
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
    
    % Get 3D dimension
    ROIlabel = reshape(ROIlabel ,256,256,256);
    [ii jj kk]=ind2sub(size(ROIlabel), find(ROIlabel==1));
    
    % Create profiles based on ROI direction
    if r==2 % Superior to inferior
        profile_SD=[]; profile_R1=[]; profile_CT=[]; profile_CU=[];
        for ll=min(jj):1:max(jj)
            profile_SD = [profile_SD; nanmean(SD_tempdata(find(jj==ll),:),1)];
            profile_R1 = [profile_R1; nanmean(R1_tempdata(find(jj==ll),:),1)];
            profile_CT = [profile_CT; nanmean(CT_tempdata(find(jj==ll),:),1)];
            profile_CU = [profile_CU; nanmean(CU_tempdata(find(jj==ll),:),1)];
        end
    elseif r==1 || r==4 || r==8  % Left to right
        profile_SD=[]; profile_R1=[]; profile_CT=[]; profile_CU=[];
        for ll=min(ii):1:max(ii)
            profile_SD = [profile_SD; nanmean(SD_tempdata(find(ii==ll),:),1)];
            profile_R1 = [profile_R1; nanmean(R1_tempdata(find(ii==ll),:),1)];
            profile_CT = [profile_CT; nanmean(CT_tempdata(find(ii==ll),:),1)];
            profile_CU = [profile_CU; nanmean(CU_tempdata(find(ii==ll),:),1)];
        end
    else % Anterior to posterior
        profile_SD=[]; profile_R1=[]; profile_CT=[]; profile_CU=[];
        for ll=min(kk):1:max(kk)
            profile_SD = [profile_SD; nanmean(SD_tempdata(find(kk==ll),:),1)];
            profile_R1 = [profile_R1; nanmean(R1_tempdata(find(kk==ll),:),1)];
            profile_CT = [profile_CT; nanmean(CT_tempdata(find(kk==ll),:),1)];
            profile_CU = [profile_CU; nanmean(CU_tempdata(find(kk==ll),:),1)];
        end
    end
    
    % Store normalized profiles
    ALL_profiles_SD{r} = normalize(profile_SD);
    ALL_profiles_R1{r} = normalize(profile_R1);
    ALL_profiles_CU{r} = normalize(profile_CU);
    ALL_profiles_CT{r} = normalize(profile_CT);
    
    % Create table for LME model
    L_ROI = size(ALL_profiles_SD{r},1); 
    
    tb{r} = [];
    frame1 = [];
    for sub=1:length(subj)
        frame1 = [frame1 repmat(group(sub), 1, L_ROI)];
    end
    
    frame2 = reshape(ALL_profiles_SD{r},1,L_ROI*length(subj));
    frame3 = reshape(ALL_profiles_CT{r},1,L_ROI*length(subj));
    frame4 = reshape(ALL_profiles_CU{r},1,L_ROI*length(subj));
    frame5 = reshape(ALL_profiles_R1{r},1,L_ROI*length(subj));
    frame6 = repmat(r,1,L_ROI*length(subj));
    
    tb{r} = table(frame1', frame2', frame3', frame4',frame5',frame6',...
            'VariableNames',{'Baby', 'SD','CT','CU','R1','Sulcus'});
    
    % Fit LME model for each sulcus
    lme{r} = fitlme(tb{r},'SD ~ 1 + CT + CU + R1 + (1|Baby)');
end

% Extract beta coefficients for plotting
betable = [];
frame1 = []; % betas
frame2 = []; % SEerr
frame3 = []; % parameter names (CT, CU, R1)
param = {'CT';'CU';'R1'};
frame4 = []; % ROI names

for r=1:length(roi_list)
    for betas=1:3
        frame1 = [frame1; lme{r}.Coefficients.Estimate(1+betas)];
        frame2 = [frame2; lme{r}.Coefficients.SE(1+betas)];
        frame3 = [frame3; {param{betas}}];
        frame4 = [frame4; {roi_list{r}}];
    end
end

betable = table(frame1, frame2, frame3, frame4,'VariableNames',{'EST', 'SE','Param','Sulcus'});

% Save beta values to CSV
writetable(betable, ['/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/data_csv_Figure_5c/', hemi, '_normalized_betas.csv');