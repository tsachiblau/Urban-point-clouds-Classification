function [pointLabels]= ...
    regionGrowingSegmentation(cloud,curvature,k,...
    angleThreshold,curvatureThreshold,minClusterSize,minDistance)

%%%%%
% Region growing segmentation (based on PCL implementation)
%
% Inputs:
%   cloud: point cloud given as MATLAB's pointCloud object. Contains N
%          points and includes normal vectors
%   k: scalar integer. Number of neighbors to consider for each point.
%   angleThreshold: scalar. Upper threshold on angle between normals. Given
%                   in radians.
%   curvatureThreshold: scalar. Upper threshold on curvature.
%   minClusterSize: scalar integer. Minimal number of points in a cluster.
%
% Outputs:
%   pointLabels: N x 1 array. Each element indicates to which cluster the
%                corresponding point belongs. Points not belonging to any
%                cluster are marked with nan.

%% Initialization
numPts = cloud.Count;
pointLabels = nan(numPts,1);
pointVisited = false(numPts,1);

% Sort curvature values and rearrange cloud
[curvatureSorted,sortCurvatureIndices] = sort(curvature,'ascend');
xyz = cloud.Location(sortCurvatureIndices,:);
normals = cloud.Normal(sortCurvatureIndices,:);

%% Find nearest neighbors of all points
[IDX_neighbors_all,distances_neighbors_all] = knnsearch(xyz,xyz,'k',k+1);

% remove self
IDX_neighbors_all(:,1) = [];
distances_neighbors_all(:,1) = [];
%% Region growing
clusterCount = 0;
numVisitedPoints = 0;

while numVisitedPoints < numPts
    
    %% Initialize new region with unvisited point with lowest curvature
    seeds = find(pointVisited==false,1);
    clusterCount = clusterCount + 1;

    % mark seed point as visited
    pointVisited(seeds) = true;
    numVisitedPoints = numVisitedPoints + 1;

    while ~isempty(seeds)
        
        %% Get current seed point coordinates, normal vector, and curvature
        IDX_currentSeed = seeds(1);        
        normal_currentSeed = normals(IDX_currentSeed,:);
      
        %% Get neighbors of current seed point
        IDX_neighbors = IDX_neighbors_all(IDX_currentSeed,:)';
        tmpDistances=distances_neighbors_all(IDX_currentSeed,:)';
        smallTmpDistances=tmpDistances<minDistance;
        IDX_neighbors=IDX_neighbors(smallTmpDistances);
        % remove visited points from neighbors list
        nbrsVisited = pointVisited(IDX_neighbors);
        IDX_neighbors(nbrsVisited) = [];
        
        % get neighbor normal vectors and curvature values
        normals_neighbors = normals(IDX_neighbors,:);
        curvature_neighbors = curvatureSorted(IDX_neighbors,:);
        
        %% Compute angle between normal vectors
        anglesBetweenNormals = acos(normals_neighbors*normal_currentSeed');
        
        %% Add points to current cluster and mark them as visited
        validNormals = or(anglesBetweenNormals < angleThreshold,...
            pi - anglesBetweenNormals < angleThreshold);        
        curvature_neighbors = curvature_neighbors(validNormals);

        IDX_addToCluster = IDX_neighbors(validNormals);
        
        pointVisited(IDX_addToCluster) = true;
        numVisitedPoints = numVisitedPoints + numel(IDX_addToCluster);

        pointLabels(IDX_addToCluster) = clusterCount;
        
        %% Add low curvature points to seed list 
        if ~isempty(IDX_addToCluster)
            validCurvature = curvature_neighbors  < curvatureThreshold;
            
            IDX_addToSeedList = ...
                IDX_addToCluster(validCurvature);
            
            seeds = [seeds; IDX_addToSeedList];
        end

        %% Remove current seed from seed list
        seeds(1) = [];
        pointLabels(IDX_currentSeed) = clusterCount;
                
    end
    
end

%% Mark points in clusters that are too small as unlabeled (nan)
numClusters = max(pointLabels);
[clusterSizes,~] = histcounts(pointLabels,(1:(numClusters+1)) - 0.5);
% histogram(pointLabels)
smallClusters_clusterNumber = find(clusterSizes < minClusterSize);

% change labeles of points belonging to small clusters to nan
Lia = ismember(pointLabels,smallClusters_clusterNumber);
pointLabels(Lia) = nan;

% rearrange labels according to curvature sorting (maintian original
% cloud's point order)
pointLabels(sortCurvatureIndices) = pointLabels;

return
