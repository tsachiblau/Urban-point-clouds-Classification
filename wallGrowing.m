function [ newPointLable,newClustersTagAndPCAndCenters ] =...
    wallGrowing( xyz,pointLable,angleThreshold,...
    distanceBetweenClusters,NumOfNeighbors,frameThickness,...
    frameExtraLength,maxDistanceForNanPoint,...
    minArea,distanceBetweenPoints,wallRecVoxelSize,dialateRadius)

    %closeClustersList=connectedClusters(xyz,pointLable,connectedDistance,connectedK,distanceBetweenPoints);
    tmpIdx=~isnan(pointLable);
    idxOfWalls=find(tmpIdx);
    
    numOfPixels=floor((1/wallRecVoxelSize)*dialateRadius);
    closeClustersList=seperatePointsToDiffrentWalls(xyz,xyz(tmpIdx,:),wallRecVoxelSize,numOfPixels);
    [~,numOfConnectedClusters]=size(closeClustersList);
    
    newPointLable=pointLable;
    
    
    
%%    
    

   
    connectedClustersList={};
    %parfor
    %get connected components
    parfor i=1:numOfConnectedClusters
        
    clusterIdx=closeClustersList{i};         
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %print group before merge
%         figure;
%         tmpLable=newPointLable(clusterIdx);
%         tmpCloud=pointCloud(xyz(clusterIdx,:));
%         tmpColored = colorPointCloud(tmpCloud,tmpLable);
%         pcshow(tmpColored);
%         xlabel('x');ylabel('y');zlabel('z');
%         title('merged walls');
%         daspect([1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        connectedClustersList{i} =connectWalls( xyz(clusterIdx,:),newPointLable(clusterIdx),...
            angleThreshold,distanceBetweenClusters,NumOfNeighbors,distanceBetweenPoints);
        

        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %print group after merge
%     figure;
%     tmpGroups=connectedClustersList{i};
%     tmpGroups=tmpGroups{1};
%     tmpIdx=find(ismember(newPointLable,tmpGroups));
%     newPointLable(tmpIdx)=20;
%     
%     tmpLable=newPointLable(clusterIdx);
%     tmpCloud=pointCloud(xyz(clusterIdx,:));
%     tmpColored = colorPointCloud(tmpCloud,tmpLable);
%     pcshow(tmpColored);
%     
%     xlabel('x');ylabel('y');zlabel('z');
%     title('merged walls');
%     daspect([1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
%% change tags


    [~,sizeOfconnectedClustersList]=size(connectedClustersList);
    for connectedClustersNum=1:sizeOfconnectedClustersList
        groupsToUnite=connectedClustersList{connectedClustersNum};
        [~,sizeOfgroupsToUnite]=size(groupsToUnite);
        for GroupsToUniteNum=1:sizeOfgroupsToUnite
            %num of diffrent groups to join
            [~,numOfNewGroups]=size(groupsToUnite{GroupsToUniteNum});

            %loop and change tags and for each cell unite the groups
            for i=1:numOfNewGroups

                %get the array of groups to join
                arrayOfUniteGroups=groupsToUnite{GroupsToUniteNum}{i};

                %get the size of the group to join
                [sizeOfGroups,~]=size(arrayOfUniteGroups);

                %if theere is no need to join
                if sizeOfGroups==1
                    continue
                end

                %
                for groupNum=2:sizeOfGroups
                    %get the idx that we need to change
                    indexToChange=newPointLable==arrayOfUniteGroups(groupNum);
                    newPointLable(indexToChange)=arrayOfUniteGroups(1);
                end  

            end
        end
    end
%% remove small walls
    newPointLable = AreaCalcAndRemove(...
    xyz,newPointLable ,minArea );

%% save normals and..
    clustersTagAndPCAndCenters = PCAOfClustersAndCenter(xyz,newPointLable); 
    newClustersTagAndPCAndCenters = clustersTagAndPCAndCenters;

%% add nan to walls.
    
    %get all the clusters 
    allClusters=clustersTagAndPCAndCenters{1};
    
    %get all the egenvectors
    clusterEgenVectors=clustersTagAndPCAndCenters{2};
    
    %get size of allClusters
    [sizeOfAllClusters,~]=size(allClusters);
    
    for i=1:sizeOfAllClusters
        currentClusterNum=allClusters(i);
        %get all the nan points
        nanIdx=find(isnan(newPointLable));
        
        %get the cluster points
        currentClusterPointsIdx=find(newPointLable==currentClusterNum);
        
        %project all cluster point on new axis 
        clusterProjection=xyz(currentClusterPointsIdx,:)*clusterEgenVectors(:,:,i);
        
        %up direction 
        %up=[0 0 1];
        
        %upProjection=up*clusterEgenVectors(:,:,i);
        %project all nan point on new axis 
        nanPointsProjection=xyz(nanIdx,:)*clusterEgenVectors(:,:,i);
        
        % get min max
        minOfCluster=min(clusterProjection);
        maxOfCluster=max(clusterProjection);
        
        %get the points in between
        boolGoodIdx=nanPointsProjection(:,1)<maxOfCluster(1)+frameExtraLength & nanPointsProjection(:,2)<maxOfCluster(2)+frameExtraLength;
        boolGoodIdx=boolGoodIdx & nanPointsProjection(:,1)>minOfCluster(1)-frameExtraLength & nanPointsProjection(:,2)>minOfCluster(2)-frameExtraLength;
        boolGoodIdx=boolGoodIdx & nanPointsProjection(:,3)>minOfCluster(3)-frameThickness & nanPointsProjection(:,3)<maxOfCluster(3)+frameThickness;
        goodIdx=find(boolGoodIdx);
        
        %we need to search for all goodIdx points that close enogh to the
        %cluster
        
        %check
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%         figure;
%         idxOfXYZ=nanIdx(goodIdx);
%         pcshow(xyz(idxOfXYZ,:),[1 0 0]);
%         xlabel('x');ylabel('y');zlabel('z');
%         title('merged walls');
%         daspect([1 1 1]);
%         hold on;
%         pcshow(xyz(currentClusterPointsIdx,:),[0 1 0]);

        %close all;
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        tmp_Idx=[];
        [numberOfnan,~]=size(goodIdx);
        Mdl = KDTreeSearcher(xyz(currentClusterPointsIdx,:));
        for goodIdxNumber=1:numberOfnan
            
            %nanIdx is points in xyz with nan lable
            %goodIdx is nan points that is close to the cluster
            idxOfXYZ=nanIdx(goodIdx(goodIdxNumber));
            [~ , D] = knnsearch(Mdl,xyz(idxOfXYZ,:),'k',1);
            [minDist,~] = min(D);
            if minDist<maxDistanceForNanPoint
                newPointLable(idxOfXYZ)=currentClusterNum;
                %needed for check
                tmp_Idx=[tmp_Idx,idxOfXYZ];
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%         figure;
%         idxOfXYZ=nanIdx(goodIdx(ismember(goodIdx,tmp_Idx) ));
%         pcshow(xyz(idxOfXYZ,:),[1 0 0]);
%         xlabel('x');ylabel('y');zlabel('z');
%         title('merged walls');
%         daspect([1 1 1]);
%         hold on;
%         pcshow(xyz(currentClusterPointsIdx,:),[0 1 0]);
% 
%         close all;
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    end
    
    
end





