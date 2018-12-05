function [ PlannarPointLabels ] = leaveOnlyWallsTags( xyz ,...
    PlannarPointLabels, clustersTagAndPC, angelTH)

    xyPlaneNormal = [0 0 1]';

    [numOfClusters,~] = size(clustersTagAndPC{1});
    eigenVectorsMat=clustersTagAndPC{2};
    allTheClusters=clustersTagAndPC{1};
    for i = 1:numOfClusters
        currCluserNormal = eigenVectorsMat(:,3,i);
        clusterTag=allTheClusters(i);
        anglesBetweenNormals = acos(currCluserNormal'*xyPlaneNormal);

        validNormals = abs((pi/2)-anglesBetweenNormals) < angelTH; 
        if ~validNormals
            PlannarPointLabels(PlannarPointLabels == clusterTag) = nan;
            % do we need to nan the clustersTagAndPC?
        end

    end
end
