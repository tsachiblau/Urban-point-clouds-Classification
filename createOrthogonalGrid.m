function [XY_indices,gridSize,xyMin,xyMax] = createOrthogonalGrid(cloud,pixelSize)

% Assumption:  gravitiy is aligned with negative z-axis direction
xy = cloud.Location(:,1:2);

xyMax = [cloud.XLimits(2) cloud.YLimits(2)];
xyMin = [cloud.XLimits(1) cloud.YLimits(1)];

gridSize = floor((xyMax-xyMin)/pixelSize) + 1;

xy_shifted = bsxfun(@minus,xy,xyMin);
XY_indices = floor(xy_shifted/pixelSize) + 1;

% z = cloud.Location(:,3);
%
% heightMap = ...
%    accumarray([XY_indices(:,2) XY_indices(:,1)],z,fliplr(gridSize),@min,fillVal;

return
