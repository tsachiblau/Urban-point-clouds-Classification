function cloud_colored = colorPointCloud(cloud,labels)

%% Inputs
% cloud - Point cloud with N points given as MATLAB's pointCloud object
% labels - N x 1 array with point labels. For unlabeled points use NaN.

%% Outputs
% cloud_colored - Point cloud with N points given as MATLAB's pointCloud 
%                 object. Each point is colored according to its label.
%                 Unlabeled points are colored in black.

%% Find unique labels
labeled = ~isnan(labels);
[labels_unique,~,label_indices] = unique(labels);

numLabels = sum(~isnan(labels_unique)); % not including unlabeled points

% each label is given an index
label_indices(label_indices>numLabels) = [];

%% Get colormap and sample it (colormap name is 'hsv')
color_map = hsv;
color_map_size = size(color_map,1);

% uniformly sample colormap
step = max([floor(color_map_size/numLabels) 1]);
sampleInds = 1:step:((numLabels-1)*step + 1);

sampleInds = 1 + mod(sampleInds,color_map_size);

color_map_sampled = color_map(sampleInds,:);

%% Assign a color to each point according to its label 
% (unlabeled points are colored in black)

% Initialize output color array
N = cloud.Count;
labelsColors = zeros(N,3);

% color points
labelsColors(labeled,:) = ...
    color_map_sampled(label_indices,:);
labelsColors = uint8(255*labelsColors);

%% Set output cloud with colors
cloud_colored = cloud;
cloud_colored.Color = labelsColors;

return
