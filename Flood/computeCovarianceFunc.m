function covarianceMatrix = computeCovarianceFunc(xyz,neighborhoodSize)
xyz = squeeze(xyz);
% pcshow(xyz)
covarianceMatrix = (xyz'*xyz)/neighborhoodSize;
return