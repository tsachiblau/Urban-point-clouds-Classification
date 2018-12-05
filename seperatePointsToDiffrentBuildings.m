function [seperateBuildings,groupsLables] = seperatePointsToDiffrentBuildings(xyz...
                                            ,xyzBuildingsAndRoof,mergedWallsPointLabels,...
                                            mergedRoofPointLabels,voxelSize,numOfPixels)
    groupsLables=zeros(length(mergedWallsPointLabels),1);
    seperateBuildings={};
    xy=xyz(:,1:2);
    
    %show 2D pictur of all walls and buildings
    smashedxyzBuildingsAndRoof=xyzBuildingsAndRoof;
    smashedxyzBuildingsAndRoof(:,3)=0;
    
    %show pic befor dialate
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
    
%     %show pic before dialate
%    	figure;
%     imagesc(binaryPic);

    
    filled_binary_pic=imfill(binaryPic);
    [lables,numOfGroups]=bwlabel(filled_binary_pic,8);

    binaryPicDilate = pumpShapes(lables,numOfPixels);
       
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %put number on the clusters 
%      figure;
%      %binaryPicDilate=binaryPic;
%      imagesc(binaryPicDilate);
%      hold on;
%     [row,column]=size(binaryPicDilate);
%     numberList=[];
%     for i=1:row
%         for j=1:column
%             if ~ismember(binaryPicDilate(i,j),numberList)
%                 text(j,i,num2str(binaryPicDilate(i,j)) ,'FontSize',10); 
%                 numberList=[numberList,binaryPicDilate(i,j)];
%             end
%         end
%      end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    tmpIterator=1;
    %sepeperate all the points to groups
    for i=1:numOfGroups
        %get the idx of the place with lable i
        tmpIdxInBinaryImage=find(binaryPicDilate==i);
        
        %how many points there are in the cluster?
        tmpIdxInClusteresXYZ=find(ismember(ind_BuildingsAndRoofInds,tmpIdxInBinaryImage));
        tmpSum=size(tmpIdxInClusteresXYZ,1);
        if tmpSum<10
            continue;
        end
        
        %find the idx the points that fit the specific cell
        tmpIdxInXYZ=find(ismember(ind_xy,tmpIdxInBinaryImage));
        tmpSum=size(tmpIdxInXYZ,1);
        if tmpSum<50
            continue;
        end
        seperateBuildings{tmpIterator}=tmpIdxInXYZ;
        
        %need to change tags to i
        tmpNeedToChangeIdx=~isnan(mergedWallsPointLabels(tmpIdxInXYZ));
        tmpNeedToChangeIdx=or(tmpNeedToChangeIdx,~isnan(mergedRoofPointLabels(tmpIdxInXYZ)));
        
        %idx to change
        tmpIdx=tmpIdxInXYZ(tmpNeedToChangeIdx);
        
        groupsLables(tmpIdx)=tmpIterator;
        tmpIterator=tmpIterator+1;
    end
    
    
    
    
%check    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tmpPointsCloud=pointCloud(xyz);    
%     for i=1:numOfGroups
%         tmpIdx=seperateBuildings{i};
%         tmpLables=zeros(length(xyz),1);
%         tmpLables(tmpIdx)=groupsLables(tmpIdx);
%         tmpLables(tmpIdx(groupsLables(tmpIdx)==0))=numOfGroups+1;
%         cloudColor= colorPointCloud( tmpPointsCloud,tmpLables);
%         figure;
%         pcshow(cloudColor);
%         xlabel('x');ylabel('y');zlabel('z');
%         title('roof and walls');
%         daspect([1 1 1]);
%         close all;
%      end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
    
end