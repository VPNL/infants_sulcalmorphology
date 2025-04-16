%% Figure 5A

morph_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/';
addpath(morph_path)

subj = {'bb28_mri12_mask'}; 

roi_list = {'calcarine','pos','insula','central','cos','sts','sfs','ips','loc','ifs','ots_lat', 'its'};

hemi = {'lh','rh'};
range1= 0;
range2= 9.5;
flag = 0;

%%

FSdir= ('/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains');
subjdir = ('/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains');

for h = 1:length(hemi)
    for l = 6 %1:length(roi_list)
        for i=1:length(subj)
            cd(FSdir)
            labelfile = readFileNifti(fullfile(FSdir, 'fsAverage', 'label','sarah_coarse_labels',[hemi{h} '.' roi_list{l} '.coarse.st.label.nii.gz']));
            cd(fullfile(subjdir, subj{i}, 'surf/'));
            measure = readFileNifti([subj{i} '_sulc_fsaverage_' hemi{h} '.nii.gz']);
            
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
            
            % Change the values to visualize each sulcus
            view([-83.6356  -7.6107])
            axis off;
            box off;
            grid off;
            clear r v c NEW COS sulc
            
            set(gcf, 'Position', [1297 320 622 501]);
            filename = fullfile(morph_path, 'figures_updated', [hemi{h} '_' subj{i} '_sulcal_depth_3d_scatter_' roi_list{l} '.eps']);
            saveas(gcf,filename,'epsc');
        end
    end
end