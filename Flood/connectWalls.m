function [groupsToUnite]=connectWalls( xyz,pointLable,angleThreshold,...
        distanceBetweenClusters,NumOfNeighbors,distanceBetweenPoints)

    
    
    %function input
    %xyz - cluster number
    %pointLable- cell array each cell have cluster number and eagenvector matrix
    %angleThreshold - below this Angle we cluster two clusters togather
    %frameThickness - the thickness of the frame that add nans
    %frameExtraLength - the extra length of the frame that add nans
    
    %function output
    %newPointLable- return new labels
    
    %initiate 
    
    newPointLable=pointLable;

    
    
%%
    %prepare for inner loop
    clustersTagAndPCAndCenters = PCAOfClustersAndCenter(xyz,newPointLable);

    %update clusters number
    clustersNormanls= squeeze(clustersTagAndPCAndCenters{2}(:,3,:));
    normalAngles=real(acos(clustersNormanls'*clustersNormanls));
    goodNormals=or(normalAngles<angleThreshold ,pi-normalAngles<angleThreshold);
    clusterList=clustersTagAndPCAndCenters{1};
    clustersWithSameNormals=[];
    [numOfClusters,~]=size(clusterList);

    
    
    %run over the whole clusters list
    newClusteridx=1;
    for i=1:numOfClusters
        
        %the current cluster
        currentCluster=clusterList(i);
        
        %find the idx of the clusters with same normal as currentCluster
        %i is the row of the current clustr
        currentClusterFriendIdx=find(goodNormals(i,:) );
        
        %if we didnt found idx then go to next iteration
        [~,notEmpty]=size(currentClusterFriendIdx);
        if notEmpty==0
           continue 
        end
        
        %update the normal list of the current group
        clustersWithSameNormals{newClusteridx}=clusterList(currentClusterFriendIdx);
                
        %update the needToCheckClusters contain the  idx of the same clusteres
        needToCheckClustersIdx=currentClusterFriendIdx;
        
        %remove the current cluster from need to check
        idxToRemove=find(needToCheckClustersIdx==currentCluster);
        needToCheckClustersIdx(idxToRemove)=[];
        
        %update the googNormals matrix so we wont check again 
        goodNormals(i,:)=0;
        goodNormals(:,i)=0;
        
        %now we ready to add the rest of the needToCheckClustersIdx array
        [~,sizeOfneedTocheckCkustersIdx]=size(needToCheckClustersIdx);
        
        while sizeOfneedTocheckCkustersIdx~=0
            
            %update the idx of the first cluster in our list
            IdxOfCluster=needToCheckClustersIdx(1);
           

            %find the idx of the clusters with same normal as currentCluster
            %i is the row of the current clustr
            currentClusterFriendIdx=find(goodNormals(IdxOfCluster,:) );
            
            %remove ethe current idx
            idxOfCurrent=find(currentClusterFriendIdx==IdxOfCluster);
            
            %if found
            if length(idxOfCurrent)>0
                currentClusterFriendIdx(idxOfCurrent)=[];
            end
            
            %if we didnt found idx then go to next iteration
            [~,notEmpty]=size(currentClusterFriendIdx);
            if notEmpty==0
                goodNormals(IdxOfCluster,:)=0;
                goodNormals(:,IdxOfCluster)=0;
                needToCheckClustersIdx(1)=[];
                [~,sizeOfneedTocheckCkustersIdx]=size(needToCheckClustersIdx);
                continue 
            end

            %update the normal list of the current group
            clustersWithSameNormals{newClusteridx}=[clustersWithSameNormals{newClusteridx} ; clusterList(currentClusterFriendIdx)];

            %remove the current cluster from need to check
            needToCheckClustersIdx(1)=[];
            
            %add cluster idx needToCheckClustersIdx list if it dosent in it
            needToCheckClustersIdx=[needToCheckClustersIdx,currentClusterFriendIdx];
            
            %remove clones
            needToCheckClustersIdx=unique(needToCheckClustersIdx);
            
            %update the googNormals matrix so we wont check again 
            goodNormals(IdxOfCluster,:)=0;
            goodNormals(:,IdxOfCluster)=0;
            
            [~,sizeOfneedTocheckCkustersIdx]=size(needToCheckClustersIdx);
        end
        clustersWithSameNormals{newClusteridx}=unique(clustersWithSameNormals{newClusteridx});
        
        %the index of the new cluster list that we are creating
        newClusteridx=newClusteridx+1;
    end
           
    
  
    
%%
    %merge all the groups 
    %parallel loop
    
    [~,sizeOfClustersWithSameNormals]=size(clustersWithSameNormals);
    groupsToUnite={};
    %parfor iteratotOfUnite=1:sizeOfClustersWithSameNormals
    parfor iteratotOfUnite = 1:sizeOfClustersWithSameNormals
        clusterNumbers=clustersWithSameNormals{iteratotOfUnite};
        %get size
        [groupSize,~]=size(clusterNumbers);
        if groupSize<1
           continue 
        end     
        
        egenVectorMat=zeros(3,3,groupSize);
        centersOfMass=zeros(groupSize,3);
        %arrange the egen vector mat
        iterationNum=0;
        for tmpClusterNum=clusterNumbers'
            iterationNum=iterationNum+1;
            clusterIdx=find(clusterList==tmpClusterNum);
            egenVectorMat(:,:,iterationNum)=clustersTagAndPCAndCenters{2}(:,:,clusterIdx);
            centersOfMass(iterationNum,:)=clustersTagAndPCAndCenters{3}(clusterIdx,:);
        end
        
        %merge close clusters
        groupsToUnite{iteratotOfUnite}=...
            MergeTheNearByClusters(xyz,newPointLable,clusterNumbers,...
            egenVectorMat,centersOfMass,...
            distanceBetweenClusters,NumOfNeighbors,distanceBetweenPoints); 
        
        
    end

    
end





