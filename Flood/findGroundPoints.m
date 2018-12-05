function groundPointsFlags = ...
    findGroundPoints(cloud,XY_indices,minMap,groundBinary,maxDiff)

[M,N] = size(groundBinary);
linearInds = sub2ind([M,N],XY_indices(:,2),XY_indices(:,1));
z = cloud.Location(:,3);

notGround = (z - minMap(linearInds)) > maxDiff;
possiblyGround = logical(groundBinary(linearInds));

groundPointsFlags = and(possiblyGround,not(notGround));

%%
% xyz = cloud.Location;
% 
% xyzGround = xyz(groundPointsFlags,:);
% pcshow(xyzGround);
% 
% xyzNotGround = xyz(not(groundPointsFlags),:);
% pcshow(xyzNotGround);

return
