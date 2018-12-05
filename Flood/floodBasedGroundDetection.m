function groundDetectionResults = ...
    floodBasedGroundDetection(cloud,pixelSize,blockSideLength,...
    elevationAngleThresh,maxPointHeightDiff,visualize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detects points belonging to ground in a point cloud using a "flood"
%   process. It is assumed that the direction of gravity in the cloud is -z
%
% Inputs:
%   cloud - point cloud given as a MATLAB pointCloud object with N points
%   pixelSize - scalar. Side length of height-map pixels. [m]
%   blockSideLength - scalar. Max. side length of largest building in the
%                     scene. [m]
%   elevationAngleThresh - scalar. Max. slope angle of ground surface in
%                          the scene. [rad]
%   visualize - show figurs? 0=no,1=yes
%
% Outputs:
%   groundDetectionResults - struct with the following fields:
%       minMap - min height-map
%       maxMap - max height-map
%       heightMapCreationData - struct. Contains data regarding height-map
%                               creation such as grid size, grid boundaries,
%                               and grid indices of each point in the cloud.
%       groundBinary - binary image. 1's correspond to ground pixels.
%       groundPointsFlags - N x 1 logical array. If an element is true,
%                           its correponding point was detected as ground.
%       groundCloud - pointCloud object. Only contains points detected as
%                     ground.
%       timing - struct. Contains time measurements for each step of
%                the algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters height difference threshold
% maxZDiff = pixelSize*tan(elevationAngleThresh);
maxZDiff = maxPointHeightDiff;

%% Create height maps (min,max)
% profile on
% tic;
[minMap,maxMap,heightMapCreationData] = ...
    createHeightMaps(cloud,pixelSize);
differenceMap = maxMap - minMap;
% groundDetectionResults.timing.heightMapsCreation = toc;
% profile viewer

%% Display height maps
if visualize
    displayHeightMaps(minMap,maxMap,differenceMap,cloud);
end

%% Get pixel indices of flood sources
% profile on
% tic;
blockSideLength_pixels = ceil(blockSideLength/pixelSize);
[floodSource,numBlocks_M,numBlocks_N] = ...
    findFloodSources(minMap,differenceMap,...
    blockSideLength_pixels,maxZDiff);
% groundDetectionResults.timing.findSources = toc;
% profile viewer

%% Display source points
if visualize
    displaySourcePoints(minMap,cloud,floodSource,...
        numBlocks_M,numBlocks_N,pixelSize);
end

%% Use flood algorithm
% profile on
% tic;
groundMask = ...
    flood(minMap,floodSource,pixelSize,elevationAngleThresh);
groundBinary = groundMask < inf;
% groundDetectionResults.timing.flood = toc;
% profile viewer

% imshow(groundMask,[])
% imshow(groundBinary,[])

%% Display ground pixels
if visualize
    displayGroundPixels(groundBinary,groundMask,minMap,cloud);
end

%% Find ground points
% tic;
groundPointsFlags = ...
    findGroundPoints(cloud,heightMapCreationData.XY_indices,...
    minMap,groundBinary,maxZDiff);
% groundDetectionResults.timing.findGrounPoints = toc;

%% Get ground cloud
groundCloud = select(cloud,find(groundPointsFlags));

%% Display ground points
if visualize
    displayGroundPoints(cloud,groundPointsFlags);
end

%% Set outputs
groundDetectionResults.minMap = minMap;
groundDetectionResults.maxMap = maxMap;
groundDetectionResults.heightMapCreationData = heightMapCreationData;

groundDetectionResults.groundBinary = groundBinary;
groundDetectionResults.groundPointsFlags = groundPointsFlags;
groundDetectionResults.groundCloud = groundCloud;

return
