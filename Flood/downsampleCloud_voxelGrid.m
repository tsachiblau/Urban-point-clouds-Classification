
function [PCL_downsampled,downsampleIndices,full3Didx] = downsampleCloud_voxelGrid(PCL,voxelSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downsample a point cloud using a voxel grid                 %
% 															  %
% Inputs:										         	  %
%		PCL - point cloud given as a MATLAB pointCloud object %
%		voxelSize - length of cubic voxel side (scalar)		  %
%															  %
% Outputs:                          								  %
%		PCL_downsampled - downsampled point cloud given as a  %
%					      MATLAB pointCloud object			  %
%		downsampleIndices - indices of remaining points 	  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get point cloud x,y,z coordinates
xyz = PCL.Location;

%% Compute voxel grid dimensions
minXYZ = min(xyz);
maxXYZ = max(xyz);
M_XYZ = ceil((maxXYZ - minXYZ)/voxelSize);

%% Find indices of the voxel containing each point
xyzInds = floor(bsxfun(@minus,xyz,minXYZ)/voxelSize) + 1;

%% Convert indices triplets into linear indices such that each voxel is represented by a single integer instead of three
ind = sub2ind(M_XYZ,xyzInds(:,1),xyzInds(:,2),xyzInds(:,3));

%% Pseudo-randomly select a single point in each voxel
[~,downsampleIndices,full3Didx] = unique(ind);
 
PCL_downsampled = select(PCL,downsampleIndices);

return
