%% Figure 5A
% This script visualizes the 3D structure of macro- and microstructural
% parameters in the 12 sulci. 

% Example data for one participant and one sulcus over the 4 timepoints are
% located in the data/Figure_5a folder. 

% NOTE: This script assumes the repository is located at ~/Downloads/infant_sulcalmorphology

%% Setup paths relative to repository root

repo_root = '~/Downloads/infant_sulcalmorphology';

% Define data path
data_path = fullfile(repo_root, 'data', 'Figure_5a');

%% Parameters
subj = {'example_subj'}; 

roi_list = {'calcarine','pos','insula','central','cos','sts','sfs','ips','loc','ifs','ots_lat', 'its'};

hemi = {'lh','rh'};

% Range will change based on which parameter you visualize
range1= 0;
range2= 9.5;

% Select which parameter to visualize (uncomment one)
parameter = 'sulc';      % Sulcal Depth
% parameter = 'thickness'; % Thickness
% parameter = 'curv';      % Curvature
% parameter = 'R1';        % R1

%%

for h = 1 % Left hemisphere
    for l = 6 % STS
        for i=1:length(subj)
            % Get the label file (ROI mask)
            label_file = fullfile(data_path, [hemi{h} '.' roi_list{l} '.label.nii.gz']);
            labelfile = readFileNifti(label_file);
            
            measure_file = fullfile(data_path, [subj{i} '_' parameter '_fsaverage_' hemi{h} '.nii.gz']);
            measure = readFileNifti(measure_file);
            
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