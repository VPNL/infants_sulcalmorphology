%% Figure 5A

%  NIfTI files for this example subject ('bb28_mri12_mask') are located in the `data` folder

morph_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/';
addpath(morph_path)

subj = {'bb28_mri12_mask'}; 

roi_list = {'calcarine','pos','insula','central','cos','sts','sfs','ips','loc','ifs','ots_lat', 'its'};

hemi = {'lh','rh'};

% Range will change based on which parameter you visualize
range1= 0;
range2= 9.5;

flag = 0;

%%

FSdir= ('/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains');
subjdir = ('/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains');

for h = 1 % Left hemisphere
    for l = 6 % STS
        for i=1:length(subj)
            cd(FSdir)
            labelfile = readFileNifti(fullfile(FSdir, 'fsAverage', 'label','sarah_coarse_labels',[hemi{h} '.' roi_list{l} '.coarse.st.label.nii.gz']));
            cd(fullfile(subjdir, subj{i}, 'surf/'));

            %% Choose which parameter to visualize by uncommenting one line

            % ---------- Sulcal Depth ----------
            % measure = readFileNifti([subj{i} '_sulc_fsaverage_' hemi{h} '.nii.gz']);

            % ---------- Thickness ----------
            % measure = readFileNifti([subj{i} '_THICK_fsaverage_' hemi{h} '.nii.gz']);

            % ---------- Curvature ----------
            % measure = readFileNifti([subj{i} '_CURV_fsaverage_' hemi{h} '.nii.gz']);

            % ---------- R1 ----------
            % whichmonth_path = extractBetween(subj{i}, '_', '_');
            % whichbaby = extractBefore(subj{i}, '_');
            % cd(['/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/' whichbaby '/' whichmonth_path '/preprocessed_aligned/qmri/'])
            % measure = readFileNifti([hemi{h} '.SEIR_R1_aligned_1mm_fsaverage_gray.nii.gz']);

            % Create masked volume
            NEW = zeros(256, 256, 256);
            NEW(labelfile.data==1) = measure.data(labelfile.data==1);

            max(max(max(NEW)))
            measure.data =NEW;
            measure.fname = 'temp.nii.gz';
            writeFileNifti(measure);
            NEW = readFileNifti('temp.nii.gz');
            !rm -f 'temp.nii.gz';
            
            [r, c, v] = ind2sub(size(NEW.data), find(NEW.data));
            val=NEW.data(find(NEW.data));
            normalized_val = normalize(val);
            
            figure;
            set(gcf, 'Color', [1 1 1])
            
            scatter3(r,c,v, 200, val, 'filled'), colormap(jet); colorbar;
            caxis([range1 range2]);

            % Adjust view (set manually using get(gca, 'View'))
            view([-119.0024 0.4599])
            axis off;
            box off;
            grid off;
            clear r v c NEW COS sulc
            
            set(gcf, 'Position', [1297 320 622 501]);
        end
    end
end