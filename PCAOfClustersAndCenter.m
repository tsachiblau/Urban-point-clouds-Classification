function [clustersPCA] = PCAOfClustersAndCenter( xyz,pointLable)
    
    %newPointLable=pointLable;
    %arrange the Labeles
    %get the nan indexes to remove in the future
    isNaNIdx=isnan(pointLable);
    notNanClusters = pointLable(~isNaNIdx);
    allTheClusters=unique(notNanClusters);
    [numOfClusters,~]=size(allTheClusters);
    
    eigenVectorsMat=zeros(3,3,numOfClusters);
    centersOfMass=zeros(numOfClusters,3);
    %get the covarience matrix of every cluster
    
    
    for i=1:numOfClusters
        
        %find the idx of the current cluster
        tmpIdx=find(pointLable==allTheClusters(i));
        
        %get the size of the cluster
        [numOfPoints,~]=size(tmpIdx);
        
        %remove the mean of the cluster
        tmpCenterOfMass=mean(xyz(tmpIdx,:));
        xyzCentered=xyz(tmpIdx,:)-tmpCenterOfMass;
        
        %calc the cov
        covarianceMatrix=(xyzCentered'*xyzCentered)/numOfPoints;

        %calc the eigenVectors,eigenValues and arrange it by value
        [eigenVectors,eigenValues] = eig(covarianceMatrix);
        eigenValues = diag(eigenValues);
        
        %get the eigenVectors
        [~,minInd] = min(eigenValues);
        [~,maxInd] = max(eigenValues);
        middleInd=find(~ismember([1 2 3],[minInd ,maxInd]));
        
        
        sortedEigenVectorsMatrix=eigenVectors(:,maxInd);
        sortedEigenVectorsMatrix=[sortedEigenVectorsMatrix, eigenVectors(:,middleInd)];
        sortedEigenVectorsMatrix=[sortedEigenVectorsMatrix,eigenVectors(:,minInd)];
        
      
        %push it to cell to return in the future
        eigenVectorsMat(:,:,i)=sortedEigenVectorsMatrix;
        centersOfMass(i,:)=tmpCenterOfMass;
    end
    clustersPCA={allTheClusters,eigenVectorsMat,centersOfMass};
end

