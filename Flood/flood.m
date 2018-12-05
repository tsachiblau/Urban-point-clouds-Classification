function groundMask = ...
    flood(heightMap,source,heightMapRes,elevationAngleThresh)

%% Initialize outputs
[M,N] = size(heightMap);
groundMask = inf(M,N);

%% Initialize visit map
visitMap = false(M,N);
visitMap(heightMap==inf) = true;
% 0 - pixel not visited yet
% 1 - visited pixel

%% Create priority queue
pq = pq_create(M*N);

%% Insert source pixel to queue
numSources = size(source,1);

for i = 1:numSources
    pq_push(pq,sub2ind([M N],source(i,1),source(i,2)),0);   
    visitMap(source(i,1),source(i,2)) = true;    
end

%% Flood
pixDistancesVec = heightMapRes*[ones(1,4) sqrt(2)*ones(1,4)];
numPixelsVisited = 0;

while 1
    
    %% 1 - Pop and mark as ground   
    if  (numPixelsVisited==M*N)
        break
    end

    [idx,elevationAngle] = pq_pop(pq);
    
    if (-elevationAngle > elevationAngleThresh)
        break;
    else
        numPixelsVisited = numPixelsVisited + 1;
        groundMask(idx) = numPixelsVisited;
    end
    
    %% 2 - Compute elevation angle to 8-adjacent neighbors 
    % that weren't inserted to queue yet
    
    % get neighbors indices
    % nbrsInds = find4nbrs(idx,M,N);
    nbrsInds = find8nbrs_v2(idx,M,N);
    
    nanMask = isnan(nbrsInds);    
    nbrsInds(nanMask) = [];
    dist2nbrs = pixDistancesVec(not(nanMask));
    
    % check if neighbors were visited 
    nbrsVisitCheck = visitMap(nbrsInds);
    
    % remove visited neighbors
    nbrsInds(nbrsVisitCheck) = [];
    dist2nbrs(nbrsVisitCheck) = [];   
    
    % compute elevation to remaining neighbors
    nbrsElevation = atan2(heightMap(nbrsInds) - heightMap(idx),dist2nbrs);
    
    %% 3 - Push neighbors to queue
    if ~isempty(nbrsInds)
        for i = 1:numel(nbrsInds)
            pq_push(pq,nbrsInds(i),-nbrsElevation(i));
            
        end
        visitMap(nbrsInds) = 1;
    end
    
end

pq_delete(pq);

return
