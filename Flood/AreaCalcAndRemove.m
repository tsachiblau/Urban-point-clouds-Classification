function [ newPointLable] = AreaCalcAndRemove( xyz,pointLable ,minArea )
    clustersPCA=[];
    newPointLable=pointLable;
    %arrange the Labeles
    %get the nan indexes to remove in the future
    isNaNIdx=isnan(pointLable);
    notNanClusters = pointLable(~isNaNIdx);
    allTheClusters=unique(notNanClusters);
    [numOfClusters,~]=size(allTheClusters);
    
    %eigenVectorsMat=zeros(3,3,numOfClusters);
    %newAllTheClusters=zeros(numOfClusters,1);
    %centerOfMass=zeros(numOfClusters,3);
    clusterIdx=0;
    %get the covarience matrix of every cluster
    for i=1:numOfClusters
        
        %find the idx of the current cluster
        tmpIdx=find(pointLable==allTheClusters(i));
        
        %get the size of the cluster
        [numOfPoints,~]=size(tmpIdx);
        
        %remove the mean of the cluster
        tmpCenterOfMass=mean(xyz(tmpIdx,:));
        tmpXYZ=xyz(tmpIdx,:)-tmpCenterOfMass;
        
        %calc the cov
        covarianceMatrix=tmpXYZ'*tmpXYZ/numOfPoints;
        
        %calc the eigenVectors,eigenValues and arrange it by value
        [eigenVectors,eigenValues] = eig(covarianceMatrix);
        eigenValues = diag(eigenValues);
        
        %get the eigenVectors
        [~,minInd] = min(eigenValues);
        [~,maxInd] = max(eigenValues);
        middleInd=find(~ismember([1 2 3],[minInd ,maxInd]));
        
        
        sortedEigenMatrix=eigenVectors(:,maxInd) ;
        sortedEigenMatrix=[ eigenVectors(:,middleInd), sortedEigenMatrix];
        sortedEigenMatrix=[sortedEigenMatrix,eigenVectors(:,minInd)];
        
        %centerOfMass(i,:)=tmpCenterOfMass;
        %calc the area of the cluster
        area= areaOfShape( xyz(tmpIdx,:),sortedEigenMatrix);

        %is the area big enoght?
        if area<minArea
            newPointLable(tmpIdx)=NaN; 
        %push it to cell to return in the future
        else
            clusterIdx=clusterIdx+1;
            %eigenVectorsMat(:,:,clusterIdx)=sortedEigenMatrix;
            %newAllTheClusters(clusterIdx)=allTheClusters(i);
        end
    end
    %eigenVectorsMat(:,:,clusterIdx+1:end)=[];
    %newAllTheClusters(clusterIdx+1:end)=[];
    %clustersPCA={newAllTheClusters,eigenVectorsMat,centerOfMass};
end

