function [minMap,maxMap,heightMapCreationData] = ...
    createHeightMaps(cloud,pixelSize)

%% Define grid over the cloud
[XY_indices,gridSize,xyMin,xyMax] = createOrthogonalGrid(cloud,pixelSize);
gridSize = fliplr(gridSize);

closingRadius = 3; % [pixels], radius of structuring element used for closing

%% Get z coordinates and set fill values
z = double(cloud.Location(:,3));
fillVal = nan;

%% Create min,max height-maps
minjmax_handle = @(x) minjmax(x);

minjmaxMap = ...
    accumarray([XY_indices(:,2) XY_indices(:,1)],z,gridSize,minjmax_handle,fillVal);
minMap = real(minjmaxMap);
maxMap = imag(minjmaxMap);

% imshow(minMap,[])
% imshow(maxMap,[])

%% Identify which pixels need to be filled
dataMask = ~isnan(minMap);
% imshow(dataMask,[])

% close dataMask
se = strel('disk',closingRadius);
dataMask_closed = imclose(dataMask,se);
% imshow(dataMask_closed,[])

pixels2fill = dataMask_closed - dataMask == 1;
% imshow(pixels2fill,[])

%% Fill detected pixels (nearest neighbor interpolation)
[minMap,maxMap] = ...
    fillPixels_NN(minMap,maxMap,dataMask,pixels2fill);
% imshow(minMap,[])
% imshow(maxMap,[])

minMap(isnan(minMap)) = inf;
maxMap(isnan(maxMap)) = -inf;

%% Set outputs
heightMapCreationData.gridSize = gridSize;
heightMapCreationData.xyMin = xyMin;
heightMapCreationData.xyMax = xyMax;
heightMapCreationData.XY_indices = XY_indices;

return
