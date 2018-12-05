%%%%%%%%%%%%%%%%%
%  main script  %
%%%%%%%%%%%%%%%%% 
    close all;
    clear all; 
    addpath(genpath('.'));  


% %%  Cloud input 
    %OrigCloudFile = 'C:\Users\orenpeer\Technion\Tsachi Blau - project A\clouds\paris\paris.ply'; % Oren
    OrigCloudFile = 'cloud.ply'; %'tsachi';
    Origcloud=pcread(OrigCloudFile);
    GeoSimFlag = 0;
%%







%% Matrix input
%     %load air and ground mat
%     %addpath(genpath('C:\Users\Oren Peer\Documents\technion\project_A\clouds\GS'))
%     AirOrigCloudFile = 'RA10-15_airborneGT.mat';
%     terrestrialOrigCloudFile= 'RA10-15_terrestrialGT.mat';
%     AirGSfile = load(AirOrigCloudFile);
%     terrestrialGSfile = load(terrestrialOrigCloudFile);
%     
%     %unite air and ground clouds
%     terestrialCloud=getfield(terrestrialGSfile,'mergedTerrestrialCloud');
%     aircloud=getfield(AirGSfile,'airborneCloud');
%     GeoSimFlag=1;
%     %create new combined cloud
%     tmp=[aircloud.Location ; terestrialCloud.Location] ; 
%     Origcloud = pointCloud(tmp);
%%







%% Ground detection using flood method.
% Set parameters
addpath(genpath('.'));
visualize = 0;
pixelSize = 1; % [m], side length of height-map pixels
blockSideLength = 5; % [m], max building-side-length
elevationAngleThresh = 25*pi/180; % [rad], max ground slope angle
maxPointHeightDiff=0.2;  % [m]

% Ground detection
tic;
if GeoSimFlag == 1
    groundDetectionResultsTerestrial = floodBasedGroundDetection(terestrialCloud,...
        pixelSize,blockSideLength,...
        elevationAngleThresh,maxPointHeightDiff,visualize);
    
    groundDetectionResultsAir = floodBasedGroundDetection(aircloud,...
        pixelSize,blockSideLength,...
        elevationAngleThresh,maxPointHeightDiff,visualize);
    
    Origcloud = pointCloud([terestrialCloud.Location ; aircloud.Location]);
    groundCompleteCloudIdx = [groundDetectionResultsTerestrial.groundPointsFlags ; groundDetectionResultsAir.groundPointsFlags];
    FinalTags = zeros(Origcloud.Count,1);
    FinalTags(groundCompleteCloudIdx)=1;
    groundSeg=toc;
    display(['ground detection time = ',num2str(groundSeg),'[sec]']);
    nonGroundCompleteCloudIdx=FinalTags==0;
    goundCloud=pointCloud(Origcloud.Location(groundCompleteCloudIdx,:));
else
    groundDetectionResults = floodBasedGroundDetection(Origcloud,...
        pixelSize,blockSideLength,...
        elevationAngleThresh,maxPointHeightDiff,visualize);
    
    groundSeg=toc;
    display(['ground detection time = ',num2str(groundSeg),'[sec]']);
    % profile viewer
    
    %change the final tags vector
    FinalTags(groundDetectionResults.groundPointsFlags)=1;
    
    %mark the idx of  ground of the complete cloud
    groundCompleteCloudIdx=FinalTags==1;
    
    %mark the idx of non ground of the complete cloud
    nonGroundCompleteCloudIdx=FinalTags==0;
    %show groud points
    goundCloud=pointCloud(Origcloud.Location(groundCompleteCloudIdx,:));
end
cloud_colored = colorPointCloud(Origcloud,FinalTags);
figure;
pcshow(cloud_colored);
xlabel('x');ylabel('y');zlabel('z');
title('ground points');
daspect([1 1 1]);

%%






%% DOWNSAMPLE - down sampling the data using voxel grid based sampling.
    % section's inputs:
    % 1. VoxleGrid Size.
    % 2. original pointcloud path.

    % Section's outputs:
    % 1. a cloud object of downsampled cloud.
    % 2. indices vector that matches the points from the DS cloud to their
    %    indices in the original data.

    %non ground cloud
    nonGroundCloud=pointCloud(Origcloud.Location(nonGroundCompleteCloudIdx,:));

    tic
    voxleSize = 0.2;
    [DownSampledNonGroundCloud, DownSampledNonGroundCloudIdx,full3Didx] = downsampleCloud_voxelGrid(nonGroundCloud, voxleSize);
    downsample=toc;
    display(['downsample non Groundpoitns time = ',num2str(downsample),'[sec]']);

    %show
    figure;
    pcshow(DownSampledNonGroundCloud);
    xlabel('x');ylabel('y');zlabel('z');
    title('all points');
    daspect([1 1 1]);
%%





%% detection of planar surfaces from the non-ground downsampled cloud.
    % section's inputs:
    % 1. nonGroundCloud - non ground cloud (with normals).
    % 2. nonGroundCloudCurvatures

    % Section's outputs:  
    % 1. non ground cloud with planner surface tags (other points are tags with NaN).

    %find normals
    %areaOfcluster=4;
    %k_normalEstimation = ceil(areaOfcluster/(voxleSize*voxleSize));
    k_normalEstimation=20;
    r= 2;
    minNumOfNbrs=6;
    maxCurv = 0.2;
    viewpoint = [0 0 0];
    tic
    [normals,DScurvature] = normalEstimation_range_nlcf(DownSampledNonGroundCloud,r, minNumOfNbrs, maxCurv,viewpoint);
    normalEstimationTime = toc;
    DownSampledNonGroundCloud.Normal = normals;
    display(['Normal estimation time = ',num2str(normalEstimationTime),'[sec]']);
%%
    %A is min size of wall
    A=0.6;
    %set the values of region growing
    k_neighbors =20;
    angleThreshold = 1*pi/180;
    curvatureThreshold = 1/20;
    minClusterSize = floor(A/(voxleSize^2));
    maxDistance=0.3; %min distance between points that we want to add to region
    tic
    
    %do region growing
    PlannarPointLabels = ...
        regionGrowingSegmentation(DownSampledNonGroundCloud,DScurvature,k_neighbors,...
        angleThreshold,curvatureThreshold,minClusterSize,maxDistance);

    %get the PCA of each cluster
    clustersTagAndPCAndCenters = PCAOfClustersAndCenter(DownSampledNonGroundCloud.Location,PlannarPointLabels);

    PlannarRegionGrowingTime = toc;
    % profile viewer

    display(['detection of planar surfaces time = ',num2str(PlannarRegionGrowingTime),'[sec]']);

    % Display segmentation result
    Plannar_cloud_colored = colorPointCloud(DownSampledNonGroundCloud,PlannarPointLabels);
    figure;
    pcshow(Plannar_cloud_colored);
    xlabel('x');ylabel('y');zlabel('z');
    title('Planar Segmentation result');
    daspect([1 1 1]);
%%







%% Walls Req
    %[~, numOfClusters] = size(clustersTagAndPC);
    angelTH = 5*pi/180;
    wallsPointLabels = leaveOnlyWallsTags(DownSampledNonGroundCloud.Location, PlannarPointLabels, clustersTagAndPCAndCenters,angelTH);
    Plannar_cloud_colored = colorPointCloud(DownSampledNonGroundCloud,wallsPointLabels);

    figure;
    pcshow(Plannar_cloud_colored);

    xlabel('x');ylabel('y');zlabel('z');
    title('Walls Segmentation result');
    daspect([1 1 1]);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %put number on the clusters 
    %     hold on;
    %     tmpClustersTagAndPCAndCenters = PCAOfClustersAndCenter(nonGroundCloud.Location,wallsPointLabels);
    %     tmp_clusterNum=tmpClustersTagAndPCAndCenters{1};
    %     tmp_centers=tmpClustersTagAndPCAndCenters{3};
    %     [tmp_size,~]=size(tmpClustersTagAndPCAndCenters{1});
    %     for i=1:tmp_size
    %         text(tmp_centers(i,1),tmp_centers(i,2),tmp_centers(i,3),num2str(tmp_clusterNum(i)),'FontSize',6);
    %     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%










%% walls cluster merging 
    angleThreshold=(10/180)*pi;
    frameThickness=0.1;
    frameExtraLength=0.1;
    minArea=10;
    distanceBetweenClusters=1;
    NumOfNeighbors=30;
    tic
    %initiate connectedClusters
    
    wallDialateRadius=1.5;
    maxDistanceForNanPoint=0.5;
    distanceBetweenPoints=voxleSize/2;
    wallRecVoxelSize=0.2;
    [mergedWallsPointLabels,ClustersTagAndPCAndCenters]=wallGrowing( DownSampledNonGroundCloud.Location,...
        wallsPointLabels,angleThreshold,...
        distanceBetweenClusters,NumOfNeighbors,...
        frameThickness,frameExtraLength,...
        maxDistanceForNanPoint,minArea,...
        distanceBetweenPoints,wallRecVoxelSize,wallDialateRadius);
    wallGrowingTime=toc;
    display(['Wall Growing estimation time = ',num2str(wallGrowingTime),'[sec]']);
       
    mergedWallsColored= colorPointCloud( DownSampledNonGroundCloud,mergedWallsPointLabels );
    
    figure;
    pcshow(mergedWallsColored);
    xlabel('x');ylabel('y');zlabel('z');
    title('merged walls');
    daspect([1 1 1]);
%%





%%
    %get roof
    DownSampledNonGroundNanIdx=(isnan(mergedWallsPointLabels));
    DownSampledNonGroundNanCloud=pointCloud( DownSampledNonGroundCloud.Location(DownSampledNonGroundNanIdx,:)  );
    DownSampledNonGroundNanCloud.Normal=normals(DownSampledNonGroundNanIdx,:);
    DownSampledNonGroundNanDScurvature=DScurvature(DownSampledNonGroundNanIdx);
        %A is min size of wall
    roofClusterSize=0.6;

    %set the values of region growing
    k_neighbors =15;
    angleThreshold = 1*pi/180;
    curvatureThreshold = 1/20;
    minClusterSize = floor(roofClusterSize/(voxleSize^2));
    maxDistance=0.5;
    tic
    
    %do region growing
    roofPointLabels = ...
        regionGrowingSegmentation(DownSampledNonGroundNanCloud,...
            DownSampledNonGroundNanDScurvature,k_neighbors,...
            angleThreshold,curvatureThreshold,minClusterSize,maxDistance);
	
    roofGrowingTime=toc;
    display(['roof Growing estimation time = ',num2str(roofGrowingTime),'[sec]']);
        
    %show roof segmentation
    roofClusersColored= colorPointCloud( DownSampledNonGroundNanCloud,roofPointLabels );
    figure;
    pcshow(roofClusersColored);
    xlabel('x');ylabel('y');zlabel('z');
    title('roof region growing');
    daspect([1 1 1]);    
        
%%






%%
    %roof growing
    angleThreshold=(180/180)*pi;
    frameThickness=0.1;
    frameExtraLength=0.1;
    minArea=10;
    distanceBetweenClusters=0.5;
    NumOfNeighbors=30;
    
    %initiate connectedClusters
       
    maxDistanceForNanPoint=0.5;
    roofDialateRadius=0.5;
    distanceBetweenPoints=voxleSize/2;
    roofRecVoxelSize=0.2;
    tic
    
    [mergedRoofPointLabels,roofClustersTagAndPCAndCenters]=wallGrowing( DownSampledNonGroundNanCloud.Location,...
        roofPointLabels,angleThreshold,...
        distanceBetweenClusters,NumOfNeighbors,...
        frameThickness,frameExtraLength,maxDistanceForNanPoint,...
        minArea,distanceBetweenPoints,...
        roofRecVoxelSize,roofDialateRadius);
    
    roofGrowingTime=toc;
    display(['roof Growing estimation time = ',num2str(roofGrowingTime),'[sec]']);
    
    roofMergedColor= colorPointCloud( DownSampledNonGroundNanCloud, mergedRoofPointLabels);
    figure;
    pcshow(roofMergedColor);
    xlabel('x');ylabel('y');zlabel('z');
    title('roof growing');
    daspect([1 1 1]);
       
%%






%%  seperate to areas
    
    %get all the roofs and all the walls points
    %the wall first then the roofs
    wallsAndRoofXYZ=[DownSampledNonGroundCloud.Location(~DownSampledNonGroundNanIdx,:) ;
        DownSampledNonGroundNanCloud.Location(~isnan(mergedRoofPointLabels),:) ];
    
    %get the idx
    idxOfRealvector=(find(DownSampledNonGroundNanIdx));
    nonGroundWallsAndRoofIdx=[find(~DownSampledNonGroundNanIdx);
        idxOfRealvector(~isnan(mergedRoofPointLabels))  ];
    
    %change mergedRoofPointLabels to the size of mergedWallsPointLabels
    upgradedMergedRoofPointLabels=mergedWallsPointLabels;
    upgradedMergedRoofPointLabels(:)=nan;
    upgradedMergedRoofPointLabels(DownSampledNonGroundNanIdx)=mergedRoofPointLabels;
    
    
    seperateVoxelSize=0.2;
    dialateRadius=10;
    numOfPixels=floor((1/seperateVoxelSize)*dialateRadius);
    tic
    [seperateBuildings,groupsLables]=seperatePointsToDiffrentBuildings(...
        DownSampledNonGroundCloud.Location,wallsAndRoofXYZ,...
        mergedWallsPointLabels,upgradedMergedRoofPointLabels,seperateVoxelSize,numOfPixels);
    
    seperateToBuildings=toc;
    display(['seperate To Buildings time = ',num2str(seperateToBuildings),'[sec]']);
    groupsLables(find(groupsLables==0))=nan;
    
    figure;
    pcshow(wallsAndRoofXYZ);
    xlabel('x');ylabel('y');zlabel('z');
    title('roof and walls');
    daspect([1 1 1]);
%%




%%
%show results till now
    %with ground
    %cloudColor= colorPointCloud( 	pointCloud([goundCloud.Location;DownSampledNonGroundCloud.Location]),[ones(goundCloud.Count,1)*1000000 ; groupsLables] );
    cloudColor= colorPointCloud( 	pointCloud( DownSampledNonGroundCloud.Location),groupsLables );    
    figure;
    pcshow(cloudColor);
    xlabel('x');ylabel('y');zlabel('z');
    title('before adding nans');
    daspect([1 1 1]);
    
    %display nans
        

    nanIdx=isnan(groupsLables);
    nanCloud=pointCloud( DownSampledNonGroundCloud.Location(nanIdx,:));
    nan_k_normalEstimation=5;
    
    [~,nanDScurvature] = normalEstimation_knn(nanCloud,...
                            nan_k_normalEstimation,viewpoint);
   
    
%%  





%%   update curvature
    nanCloudIdx = isnan(groupsLables);
    nanCloud = pointCloud(DownSampledNonGroundCloud.Location(nanCloudIdx,:));
    rNan=0.5;
    minNumOfNbrsNan=10;
	[normalsNew,DScurvatureNew] = normalEstimation_range_nlcf(nanCloud,rNan, minNumOfNbrsNan, maxCurv,viewpoint);
    DScurvature(nanCloudIdx)=DScurvatureNew;
    normals(nanCloudIdx,:)=normalsNew;
    
    xyz=DownSampledNonGroundCloud.Location;
    
    %print normals
    
%     figure;
%     pcshow(xyz);
%     hold on;
    %quiver3(xyz(:,1),xyz(:,2),xyz(:,3),normalsNew(:,1),normalsNew(:,2),normalsNew(:,3),'m');
    %xlabel('x');ylabel('y');zlabel('z');
    %title('Estimated normal vectors');
    %colormap jet
    %daspect([1 1 1]);
    
    %print curvature
    
%     figure;
%     pcshow(xyz,DScurvature);
%     xlabel('x');ylabel('y');zlabel('z');
%     c = colorbar;
%     c.Label.String = 'Curvature';
%     title('Estimated curvature');
%     colormap jet
%     caxis([0 1/3]);
%     daspect([1 1 1]);
    
%%








%% add nans
    nanAddedLables=groupsLables;
    [~,numOfSeperateGroups]=size(seperateBuildings);
    tmpGroupLables={};
    
    numOfNighbors=5;
    numOfTopNighbors=5;
    
    tic
    %parfor
    parfor i=1:numOfSeperateGroups
        %i=i
        tmpGroupIdx=seperateBuildings{i};
        tmpPoints=DownSampledNonGroundCloud.Location(tmpGroupIdx,:);
        [~,I]=sort(tmpPoints(:,3));
        tmpLables=groupsLables(tmpGroupIdx);
        tmpCurveture=DScurvature(tmpGroupIdx);
        newTmpLables=pointsConnectedToCluster(tmpPoints(I,:),...
                            tmpLables(I),numOfNighbors,...
                            numOfTopNighbors,tmpCurveture);
        tmpGroupLables{i}=newTmpLables(I);
    end
    
    addNans=toc;
    display(['add Nans time = ',num2str(addNans),'[sec]']);
    
    % change tags
    for i=1:numOfSeperateGroups
        tmpGroupIdx=seperateBuildings{i};
        nanAddedLables(tmpGroupIdx)=tmpGroupLables{i};
    end

%%




%% show results till now
    %with ground
    %cloudColor= colorPointCloud( 	pointCloud([goundCloud.Location;DownSampledNonGroundCloud.Location]),[ones(goundCloud.Count,1)*1000000 ; nanAddedLables] );
    
    cloudColor= colorPointCloud( 	pointCloud(DownSampledNonGroundCloud.Location),nanAddedLables )
    figure;
    pcshow(cloudColor);
    xlabel('x');ylabel('y');zlabel('z');
    title('after adding nans');
    daspect([1 1 1]);
    
%%





%% error calc to Down Sampled
 realLables = csvread('realLables.csv'); %csvread('realLables.csv');
%  realLables = realLables(:,4);
%  save('realLables.mat','realLables');
%  csvwrite('realLables.csv',realLables);

% groupsLablesMod = groupsLables;
% groupsLablesMod(~isnan(groupsLables))= 2;
% groupsLablesMod(isnan(groupsLables))= 0;
if GeoSimFlag == 1
    realLablesTerrestrial = load('RA10-15_terrestrialGT');
    realLablesTerrestrial = realLablesTerrestrial.pointClass_terrestrial;
    
    realLablesairborne = load('RA10-15_airborneGT');
    realLablesairborne = realLablesairborne.pointClass_airborne;
    
    realLables = [realLablesTerrestrial;realLablesairborne];
    groundDetectionTags = [ groundDetectionResultsAir.groundPointsFlags ; groundDetectionResultsTerestrial.groundPointsFlags];
    nonGroundIdx = find ( ~groundCompleteCloudIdx) ;
    nanAddedLablesMod = nanAddedLables;
    nanAddedLablesMod(~isnan(nanAddedLablesMod))= 2;
    nanAddedLablesMod(isnan(nanAddedLablesMod))= 0;
    
else
    groundDetectionTags = groundDetectionResults.groundPointsFlags;
    nonGroundIdx = find ( ~groundDetectionResults.groundPointsFlags);
    nanAddedLablesMod = nanAddedLables;
    nanAddedLablesMod(~isnan(nanAddedLablesMod))= 2;
    nanAddedLablesMod(isnan(nanAddedLablesMod))= 0;
end

errorCalcDS( double(groundDetectionTags) ,nanAddedLablesMod,...
           nonGroundIdx(DownSampledNonGroundCloudIdx), realLables,...
             GeoSimFlag,DownSampledNonGroundCloud, Origcloud); %%



%  groundCompleteCloudIdx, nonGroundCompleteCloudIdx ,OurDSTags,ind )
%% upsample 
    upsampledTags = upsample_voxleGrid(groundCompleteCloudIdx, nonGroundCompleteCloudIdx , ...
        nanAddedLables,full3Didx );
    cloudColor= colorPointCloud(Origcloud , upsampledTags);
    figure;
    pcshow(cloudColor);
    xlabel('x');ylabel('y');zlabel('z');
    title('after upsampling');
    daspect([1 1 1]);

%%




%% error calc Final (!!)
   %csvread('realLables.csv');
%  realLables = realLables(:,4);
%  save('realLables.mat','realLables');
%  csvwrite('realLables.csv',realLables);

if GeoSimFlag == 1
    realLablesTerrestrial = load('RA10-15_terrestrialGT');
    realLablesTerrestrial = realLablesTerrestrial.pointClass_terrestrial;
    
    realLablesairborne = load('RA10-15_airborneGT');
    realLablesairborne = realLablesairborne.pointClass_airborne;
    
    realLables = [realLablesTerrestrial;realLablesairborne];
    groundDetectionTags = [ groundDetectionResultsAir.groundPointsFlags ; groundDetectionResultsTerestrial.groundPointsFlags];
    nonGroundIdx = find ( ~groundCompleteCloudIdx) ;
    upsampledTagsMod = upsampledTags;
    upsampledTagsMod(~isnan(upsampledTagsMod))= 2;
    upsampledTagsMod(isnan(upsampledTagsMod))= 0;
    
else
    realLables = csvread('realLables.csv');
    upsampledTagsMod = upsampledTags;
    upsampledTagsMod(~isnan(upsampledTags) & upsampledTags~=1)= 2;
    upsampledTagsMod(isnan(upsampledTags))= 0;
end
%%



%% display final result!!!
   cloudColor= colorPointCloud(Origcloud , upsampledTagsMod);
    figure;
    pcshow(cloudColor);
    xlabel('x');ylabel('y');zlabel('z');
    title('after all');
    daspect([1 1 1]);
    
    
    
%%



%% error calc final!
errorCalc(upsampledTagsMod , realLables,GeoSimFlag, Origcloud); %%
