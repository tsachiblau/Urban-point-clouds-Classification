function [seperateBuildings] = seperatePointsToDiffrentWalls(xyz...
                                            ,xyzBuildingsAndRoof,voxelSize,numOfPixels)
    seperateBuildings={};
    xy=xyz(:,1:2);
    
    %show 2D pictur of all walls and buildings
%    smashedxyzBuildingsAndRoof=xyzBuildingsAndRoof;
%    smashedxyzBuildingsAndRoof(:,3)=0;
%     figure;
%     pcshow(smashedxyzBuildingsAndRoof);
%     xlabel('x');ylabel('y');zlabel('z');
%     title('roof and walls');
%     daspect([0.1 0.1 1]);
    
    % seperate to voxels
    minXY = min(xy);
    maxXY = max(xy);
    M_XY = ceil((maxXY - minXY)/voxelSize);
    xyInds = floor(bsxfun(@minus,xy,minXY)/voxelSize) + 1;
    xyBuildingsAndRoofInds=floor(bsxfun(@minus,xyzBuildingsAndRoof(:,1:2),minXY)/voxelSize) + 1;
    ind_xy = sub2ind(M_XY,xyInds(:,1),xyInds(:,2) );
    ind_BuildingsAndRoofInds = sub2ind(M_XY,xyBuildingsAndRoofInds(:,1),xyBuildingsAndRoofInds(:,2) );
    binaryPic=zeros(M_XY);
    
    
    binaryPic(ind_BuildingsAndRoofInds)=1;
    
    %show befor dialate
%     figure;
%     imagesc(binaryPic);

    
    SE = strel('disk',numOfPixels);
    binaryPicDilate = imdilate(binaryPic , SE);
    filled_binary_pic=imfill(binaryPicDilate);
    [lables,numOfGroups]=bwlabel(filled_binary_pic,4);
    
    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %put number on the clusters 
%     figure;
%     imagesc(lables);
%     hold on;
%     [row,column]=size(lables);
%     numberList=[];
%     for i=1:row
%         for j=1:column
%             if ~ismember(lables(i,j),numberList)
%                 text(j,i,num2str(lables(i,j)) ,'FontSize',10); 
%                 numberList=[numberList,lables(i,j)];
%             end
%         end
%     end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    %sepeperate all the points to groups
    
    for i=1:numOfGroups
        %get the idx of the place with lable i
        tmpIdxInBinaryImage=find(lables==i);
        
        %find the idx the points that fit the specific cell
        tmpIdxInXYZ=find(ismember(ind_xy,tmpIdxInBinaryImage));
        seperateBuildings{i}=tmpIdxInXYZ;
    end
    
    
    
    
%check    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure;
%     for i=1:numOfGroups
%         tmpIdx=seperateBuildings{i};
%         tmpPointsCloud=pointCloud(xyz(tmpIdx,:));
%         cloudColor= colorPointCloud( tmpPointsCloud,zeros(length(tmpIdx),1));
%         figure;
%         pcshow(cloudColor);
%         xlabel('x');ylabel('y');zlabel('z');
%         title('roof and walls');
%         daspect([1 1 1]);
% 
%         close all;
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
    
end