function[connectedClustersList] =connectedClusters(xyz,pointLable,distTH, k,distanceBetweenPoints)
    

    isNaNIdx=isnan(pointLable);
    xyClustersOnlyLables = pointLable(~isNaNIdx);
       
    xyClustersOnly=xyz(~isNaNIdx,:);
    xyClustersOnly(:,3)=0;
    
    %show clusters on xy plane
    xyClusterOnlyCloud=pointCloud(xyClustersOnly);
    
    %get the 
    xyClustersTagAndPCAndCenters = PCAOfClustersAndCenter(xyClustersOnly,xyClustersOnlyLables);
    
%     xyPlane= colorPointCloud( xyClusterOnlyCloud,xyClustersOnlyLables );
%     figure;
%     pcshow(xyPlane);
%     xlabel('x');ylabel('y');
%     title('xyPlane');
    
    
    
    %connect the clusters
    connectedClustersList=MergeTheNearByClusters (...
    xyClustersOnly, xyClustersOnlyLables , ...
    xyClustersTagAndPCAndCenters{1} ,...
    xyClustersTagAndPCAndCenters{2}(:,:,:) ,...
    xyClustersTagAndPCAndCenters{3}(:,:), distTH, k,distanceBetweenPoints);

% check
%     
%     [~,sizeOfconnectedClustersList]=size(connectedClustersList);
%     cloud=pointCloud(xyz(~isNaNIdx,:));
%     for i=1:sizeOfconnectedClustersList
%         group=cell2mat(connectedClustersList(i));
%         idx=ismember(pointLable(~isNaNIdx),group);
%         
%         tmpPic= colorPointCloud( cloud,idx);
%     
%         figure;
%         pcshow(tmpPic);
%         xlabel('x');ylabel('y');zlabel('z');
%         title('merged walls');
%         daspect([1 1 1]);    
%         close all;
%     end
end