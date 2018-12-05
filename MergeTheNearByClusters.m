function [changeList] = MergeTheNearByClusters ( xyz, pointLables , ...
    tags , egenVectorMat , centersOfMass , distTH, k,distanceBetweenPoints)

% get the hull points
    
    %get the num of clusters
    [numOfClusters,~]=size(tags);
    
    %initiate loop
    actualHullPoints={};
    for i = 1:numOfClusters
        
        %get the idx of the cluster 
        clusterIdx = find (pointLables == tags(i));
        
        %get the xyz of the cluster 
        clustersXYZ = xyz( clusterIdx , : );
        
        %project the point on egenVectorMat
        projectionOfClusterXYZ=clustersXYZ*egenVectorMat(:,:,i);
        
        %get the hull of the points
        convHull = convhull( double(projectionOfClusterXYZ(:,1)) ,double(projectionOfClusterXYZ(:,2)));
        %each cell: idx of hull in orig XYZ  | the points
        actualHullPoints{i} = addDummies( xyz(clusterIdx(convHull),:), distanceBetweenPoints);
    end

%% build kd-tree to centers of masses, and try to region grow on them( using clusters points) 

    % Sort Size values and rearrange cloud
    % [ cluster num sorted , index in cluster array!]
    % [SizeSorted,sortSizeIndices] = sort(clusterSizes,'ascend');
    % clustersCOM = hullArray{sortSizeIndices};

    %Find nearest neighbors of all clusters
    
    %do the search
    [centersDistanceIdx,~] = knnsearch(centersOfMass,centersOfMass,'k',k+1);
    
    % remove self COM (center of mass)
    centersDistanceIdx(:,1) = [];

%% Region growing

    % Initialization 
    
    %we will go throw k nearest neighbors or num of clusters
    numOfNeighborsForRegionGrowing=min(numOfClusters,k);
    
    %signed all the clusters as not visited yet
    clustersVisited = false(numOfClusters,1);
    newClusterCount = 0;
    changeList={};
    %while loop
    while sum(clustersVisited) < numOfClusters

        % Initialize new region 
        
        %seed contain the idx of the clusters need to visits
        seeds = find(clustersVisited==false,1);
        
        %the idx of the new cluster that will be create
        newClusterCount = newClusterCount + 1;

        % mark seed cluster as visited
        clustersVisited(seeds) = true;
        
        AddToChangeList=[tags(seeds(1))];
        
        while ~isempty(seeds)

            % Get current seed idx of the tags vector
            currentSeedIdx = seeds(1);  
            
            currentSeedHull = actualHullPoints{currentSeedIdx} ;

            % Get center nighbors of current seed cluster
            neighborsIdx = centersDistanceIdx(currentSeedIdx,:)';

            % creat vector with 1 in every cluster that has been visited
            nbrsVisited = clustersVisited(neighborsIdx);
            
            %remove the idx that had been visited
            neighborsIdx(logical(nbrsVisited)) = []; 
            %
            tmp=0;
            tmp=sum(ismember(neighborsIdx,[93]));
            %now neighborsIdx contaid the idx of the clusters that need to
            %be visited among the clusters that havent been visited
            
            % get neighbors min distsance values
            [numOfcenterClosestSeeds,nullVar]=size(neighborsIdx);
            allMinDist=zeros(numOfcenterClosestSeeds,1); 
            
            %check if neighborsIdx contain enything
            if nullVar==0
                seeds(1) = [];
                continue;
            end
          	for i=1:numOfcenterClosestSeeds

            	%get cluster idx
                tmpClusterIdx=neighborsIdx(i);

             	%get all the distances
               	[~ , D] = knnsearch(currentSeedHull, actualHullPoints{tmpClusterIdx});

              	%the min distance between two clusters
              	[minDist,~] = min(D);
              	allMinDist(i) = minDist;

            end 

            %looking for clusters with distance smaller that distTH
            validDistances = find(allMinDist < distTH);   
            
            if validDistances
                addToClusterIdx = neighborsIdx(validDistances);
                clustersVisited(addToClusterIdx) = true;
                AddToChangeList = [AddToChangeList;tags(addToClusterIdx)];
                
                % Add all clusters to seed list 
                seeds = [seeds;addToClusterIdx];
            end
            
            % Remove current seed from seed list
            seeds(1) = [];
            
        end
        
        
        if AddToChangeList
            changeList{newClusterCount} = AddToChangeList;
        end
        
    end
    
end
