%% This script generates Figures 2, 3, and 4

% Note: SP stands for sulcal span and SD stands for sulcal depth. 

clear; clc; close all;

% Uncomment the parameter you want to run:
%param = 'AREA';  
%param = 'SULC'; 
%param = 'R1'; 
%param = 'THICK'; 
param = 'CURV'; 

% Set axis limits for figures
switch param
    case 'AREA'
        yLimModel = [0 60];
        xLimModel = [0 420];
        yLimIntercept = [0 25];
        yLimSlope     = [0 10];
    case 'SULC'
        yLimModel = [0 8];      
        xLimModel = [0 420];
        yLimIntercept = [0 7];  
        yLimSlope     = [-0.5 1.5];   
    case 'R1'
        yLimModel = [0 0.8];     
        xLimModel = [0 420];
        yLimIntercept = [0 0.4]; 
        yLimSlope     = [0 0.16]; 
    case 'THICK'
        yLimModel = [0 4];       
        xLimModel = [0 420];
        yLimIntercept = [0 2.5];   
        yLimSlope     = [-0.2 0.8];  
    case 'CURV'
        yLimModel = [0 0.25];      
        xLimModel = [0 420];
        yLimIntercept = [0 0.25];  
        yLimSlope     = [-0.03 0];  
end

% Paths to change
FSdir      = '/oak/stanford/groups/kalanit/biac2/kgs/anatomy/freesurferRecon/babybrains';
morph_path = '/oak/stanford/groups/kalanit/biac2/kgs/projects/babybrains/mri/code/morphology_allparameters/';
saveDir = fullfile(morph_path, 'REVISION', 'Figure_2');
mkdir('Figure_2')
set_participant_information 

% Variables to define
hemis    = {'lh','rh'};
numHemis = length(hemis);
numROIs  = length(roi_list);

% Store model parameters
inC2    = zeros(numHemis, numROIs);
slP2    = zeros(numHemis, numROIs);
inCSE2  = zeros(numHemis, numROIs);
slPSE2  = zeros(numHemis, numROIs);
R_values    = zeros(numHemis, numROIs);
p_values    = zeros(numHemis, numROIs);
adjusted_Rsq= zeros(numHemis, numROIs);
model1      = cell(numHemis, numROIs);

%% Loop over hemispheres
for h = 1:numHemis
    thisHemi = hemis{h};   % 'lh' or 'rh'   
    csvFile = [param '_' thisHemi '.csv'];
 
    
    T = readtable(csvFile);                     % Read the csv file
    T = removevars(T, "Row");
    param_data_full = table2array(T);           % Convert table to numeric matrix
    roiNames = T.Properties.VariableNames;      % Store the ROI names from the table

    % Create a figure for the model data (with subplots for each ROI)
    figModel = figure('Name',[param ' - ' thisHemi],'Color','white');
    hold on; 
    
    % For each ROI
    for r = 1:numROIs
        
        % Extract data for this ROI
        param_data = param_data_full(:, r); 
        
        % Fit the linear mixed effects model
        tbl = table(logage', group', param_data_full(:,r), ...
            'VariableNames', {'logage','Baby','Param'});
        lme1 = fitlme(tbl, 'Param ~ logage + (1|Baby)');
        
        % Collect coefficients
        intercept   = lme1.Coefficients.Estimate(1);
        slope       = lme1.Coefficients.Estimate(2);
        interceptSE = lme1.Coefficients.SE(1);
        slopeSE     = lme1.Coefficients.SE(2);
        
        % Save to arrays
        inC2(h,r)   = intercept;
        slP2(h,r)   = slope;
        inCSE2(h,r) = interceptSE;
        slPSE2(h,r) = slopeSE;
        
        % Also store correlation
        [R, p] = corrcoef(param_data_full(:,r), logage);
        R_values(h, r) = R(1, 2);
        p_values(h, r) = p(1, 2);
        
        % Get the variance 
        adjusted_Rsq(h, r) = lme1.Rsquared.Adjusted;
        model1{h, r} = lme1;
        
        subplot(1, numROIs, r); 
        hold on;
        
        % Scatter plot
        scatter(10.^logage, param_data_full(:,r), ...
            50, logage, 'filled', ...
            'MarkerFaceColor', color(r,:), ...
            'MarkerEdgeColor', [0.7 0.7 0.7]);
        
        % Plot the fitted line
        xVals = linspace(xLimModel(1), xLimModel(2), 200); 
        yVals = intercept + slope .* log10(xVals); 
        plot(xVals, yVals, 'Color', [0.3 0.3 0.3], 'LineWidth', 3);
        
        % Axis limits
        xlim(xLimModel);
        ylim(yLimModel);
        
        % Make figure ready
        box off
        set(gca,'color','white')
        if r==1
            set(gca,'XTick',[]) 
            yticklabels("")
        else
            set(gca,'XTick',[])
            set(gca,'YColor','none') 
        end
        hold off;
        
    end 
    
    set(figModel, 'Position', [100 300 1400 400]); 

     % Save the figure
    figModelFile = fullfile(saveDir, sprintf('Model_%s_%s.png', param, thisHemi));
    saveas(figModel, figModelFile); 
    
end 

%% Plot intercept & slope per hemisphere
for h = 1:numHemis
    
    thisHemi = hemis{h};
    
    % Intercept 
    figIntercept = figure('Name',['Intercept - ' param ' - ' thisHemi],'Color','white');
    hold on;
    set(gca,'XTick',[], 'YTick',[]);
    ylim(yLimIntercept);

    b = bar(1:numROIs, inC2(h,:), 'FaceColor','flat', 'EdgeColor','none');

    for r = 1:numROIs
        b.CData(r, :) = color(r, :);
    end 

    errorbar(1:numROIs, inC2(h, :), inCSE2(h, :), 'Color', [0.5 0.5 0.5], 'LineStyle','none');
    box off;
    set(gca, 'XTick', []);
    set(gca, 'FontName', 'Arial', 'YTick', linspace(min(yLimIntercept), max(yLimIntercept), 3))
    set(figIntercept, 'Position', [200 700 400 400]); 
    hold off;

    figInterceptFile = fullfile(saveDir, sprintf('Intercept_%s_%s_bar.png', param, thisHemi));
    saveas(figIntercept, figInterceptFile); 
    
    % Slope
    figSlope = figure('Name',['Slope - ' param ' - ' thisHemi],'Color','white');
    hold on;
    ylim(yLimSlope);

    b = bar(1:numROIs, slP2(h,:), 'FaceColor','flat', 'EdgeColor','none');

    for r = 1:numROIs
        b.CData(r, :) = color(r, :);
    end 

    errorbar(1:numROIs, slP2(h, :), slPSE2(h, :), 'Color', [0.5 0.5 0.5], 'LineStyle','none');
    set(gca, 'XTick', []);
    set(gca, 'FontName', 'Arial');

    box off;
    set(figSlope, 'Position', [650 700 400 400]); 
    hold off;
    
    figSlopeFile = fullfile(saveDir, sprintf('Slope_%s_%s_bar.png', param, thisHemi));
    saveas(figSlope, figSlopeFile); 
    
end
