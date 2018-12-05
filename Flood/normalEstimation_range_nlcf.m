function [normals,curvature] = ...
    normalEstimation_range_nlcf(cloud,r,minNumNBRS,maxValidCurvature,viewpoint)

%% Find points within r-neighborhood of each point in the cloud
xyz = cloud.Location;
IDX_inRange = rangesearch(xyz,xyz,r);
numPts = size(xyz,1);

%% Compute covariance for each neighborhood
normals = nan(numPts,3);
curvature = inf(numPts,1);

for i = 1:numPts
    IDX_inRange_i = IDX_inRange{i};
    numNBRS = numel(IDX_inRange_i);
    
    if numNBRS >= minNumNBRS
        xyz_NBRS = xyz(IDX_inRange_i,:);
        xyz_query = xyz_NBRS(1,:);
        xyz_NBRS = bsxfun(@minus,xyz_NBRS,xyz_query);
        
        covMat = (xyz_NBRS'*xyz_NBRS)/numNBRS;
        [V,D] = eig(covMat);
        [eigenValues_sorted,IDX_Sort] = sort(diag(D));
        normals(i,:) = V(:,IDX_Sort(1))';
        
        curvature(i) = eigenValues_sorted(1)/sum(eigenValues_sorted);
    end
end

%% Nearest low curvature fix
tooHighCurvature = curvature > maxValidCurvature;

validInds = find(not(tooHighCurvature));
xyzValid = xyz(not(tooHighCurvature),:);
xyzInvalid = xyz(tooHighCurvature,:);

% for each invalid point, find nearest point with low curvature
IDX = knnsearch(xyzValid,xyzInvalid);

% replace invalid normal with valid normal
normals(tooHighCurvature,:) = normals(validInds(IDX),:);

%% For points with not enough neighbors, return normal and curvature of closest valid point
tooFewNbrs = curvature==inf;
validInds = find(not(tooFewNbrs));

xyzValid = xyz(not(tooFewNbrs),:);
xyzInvalid = xyz(tooFewNbrs,:);

% for each invalid point, find nearest valid point
IDX = knnsearch(xyzValid,xyzInvalid);

% replace invalid normal with valid normal
normals(tooFewNbrs,:) = normals(validInds(IDX),:);
curvature(tooFewNbrs) = curvature(validInds(IDX));

%% Flip normals towards the viewpoint
view2pts = bsxfun(@minus,xyz,viewpoint);
flipFlags = sum(view2pts.*normals,2) > 0;
normals(flipFlags,:) = -normals(flipFlags,:);

return
