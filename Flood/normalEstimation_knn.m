function [normals,curvature] = ...
    normalEstimation_knn(cloud,k,viewpoint)

%%%%%
% Compute normal vector and curvaure for each point and its k-neighborhood
%
% Inputs:
%   cloud: point cloud given as MATLAB's pointCloud object. Contains N
%          points.
%   k: scalar integer. Number of neighbors to consider for each point.
%   viewpoint: 1 x 3 array. Estimated normals are flipped,if necessary,
%              towards the viewpoint, given as (x,y,z) triplet.
%
% Outputs:
%   normals: N x 3 array. Each row is the x,y,z elements of the normal
%            vector estimated for the corresponding point in the point 
%            cloud.
%   curvature: N x 1 array. Each element is the curvature of the
%              corresponding points k-neighborhood. Curvature is defined
%              here as the ratio between the smallest eigenvalue of the
%              scatter matrix and the sum of the eigenvalues.

%% Get point cloud data
xyz = cloud.Location;
% pcshow(xyz)
numPts = cloud.Count;

%% Initialize outputs
normals = zeros(numPts,3);
curvature = zeros(numPts,1);

%% find k nearest neighbors of each point in the cloud
% create kdtree
kdTreeObj = KDTreeSearcher(xyz,'distance','euclidean');

% get k nearest neighbours + the point itself
indices_knn = knnsearch(kdTreeObj,xyz,'k',(k+1));

%% compute covariance for each neighborhood

% get all neighborhoods
xyz_neighborhoods = reshape(xyz(indices_knn(:),:),[numPts k+1 3]);

% compute mean of each neighborhood
xyz_mean = mean(xyz_neighborhoods,2);
% pcshow(xyz_mean)

% center each neighborhood to its mean
xyz_neighborhoods = bsxfun(@minus,xyz_neighborhoods,xyz_mean);


% compute covariance matrix for each neighborhood
neighborhoods_cell = mat2cell(xyz_neighborhoods,ones(numPts,1));

computeCovarianceFunc_handle = @(xyz) computeCovarianceFunc(xyz,k+1);

covarianceMatrices = ...
    cellfun(computeCovarianceFunc_handle,neighborhoods_cell,...
    'UniformOutput',false);

%% normals and curvature calculation
for i = 1:numPts
    
    % get covariance matrix
    covarianceMatrix = covarianceMatrices{i}; 
    
    % get eigenvalues and eigenvectors
    [eigenVectors,eigenValues] = eig(covarianceMatrix);
    eigenValues = diag(eigenValues);
    
    [eigenValue3,minInd] = min(eigenValues);
    
    % get normals
    normals(i,:) = eigenVectors(:,minInd)';
    
    % get curvature
    curvature(i) = eigenValue3/sum(eigenValues);
end

%% Flip normals not pointing towards the viewpoint
view2pts = bsxfun(@minus,xyz,viewpoint);
flipFlags = (sum(view2pts.*normals,2) > 0);
normals(flipFlags,:) = -normals(flipFlags,:);

return
