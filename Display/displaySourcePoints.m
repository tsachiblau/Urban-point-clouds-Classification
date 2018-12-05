function displaySourcePoints(minMap,cloud,floodSource,...
    numBlocks_M,numBlocks_N,pixelSize1)

defaultFontSize = 16;
set(0,'DefaultAxesFontSize',defaultFontSize);
set(0, 'DefaultAxesFontName','Times');

xyMin = [cloud.XLimits(1) cloud.YLimits(1)];
xyMax = [cloud.XLimits(2) cloud.YLimits(2)];
xLimDisp = [xyMax(1); xyMin(1)];
yLimDisp = [xyMax(2); xyMin(2)];

[M,N] = size(minMap);
floodSource_xy = [(xyMin(1) + (floodSource(:,2) - 0.5)*pixelSize1) ...
    (xyMin(2) + (floodSource(:,1) - 0.5)*pixelSize1)];
floodSource_z = ...
    minMap(sub2ind([M N],floodSource(:,1),floodSource(:,2)));
floodSource_xyz = [floodSource_xy floodSource_z];

%%
xRange = xLimDisp(1) - xLimDisp(2);
yRange = yLimDisp(1) - yLimDisp(2);

xVec = xLimDisp(2) + (xRange/numBlocks_N)*[1:(numBlocks_N-1)];
yVec = yLimDisp(2) + (yRange/numBlocks_M)*[1:(numBlocks_M-1)];

zMin = cloud.ZLimits(1);

%% Display source points - min. map
sourceMarkerSize = 150;

figure;
imagesc(xLimDisp,yLimDisp,rot90(minMap,2));
set(gca,'YDir','normal');
xlabel('x[m]');ylabel('y[m]');
title(['Min. height map - ',num2str(M),'x',num2str(N)]);
c = colorbar;
c.Label.String = 'z[m]';
[cmin1,cmax1] = caxis;

hold on
scatter(floodSource_xy(:,1),floodSource_xy(:,2),sourceMarkerSize,'m*');

%% Display block lines - cloud
zLift = 250;
linWidth = 4;

figure
pcshow(cloud.Location);
hold on

for i = 1:(numBlocks_N-1) % vertical lines (along y)
    % Your two points
    P1 = [xVec(i),xyMin(2),zMin+zLift];
    P2 = [xVec(i),xyMax(2),zMin+zLift];
    
    % Their vertial concatenation is what you want
    pts = [P1; P2];
    
    % Because that's what line() wants to see
    line(pts(:,1), pts(:,2), pts(:,3),'Color',[1 0 1],'LineWidth',linWidth);
end

for i = 1:(numBlocks_M-1) % horizontal lines (along x)
    
    % Your two points
    P1 = [xyMin(1),yVec(i),zMin+zLift];
    P2 = [xyMax(1),yVec(i),zMin+zLift];
    
    % Their vertial concatenation is what you want
    pts = [P1; P2];
    
    % Because that's what line() wants to see
    line(pts(:,1), pts(:,2), pts(:,3),'Color',[1 0 1],'LineWidth',linWidth);
    
end


xlabel('x[m]');ylabel('y[m]');zlabel('z[m]');
title('Block division');

c = colorbar;
c.Label.String = 'z[m]';
caxis([cmin1,cmax1]);
view(0,90);

%% Display source points - cloud
figure
pcshow(cloud.Location);
hold on

for i = 1:(numBlocks_N-1) % vertical lines (along y)
    % Your two points
    P1 = [xVec(i),xyMin(2),zMin+zLift];
    P2 = [xVec(i),xyMax(2),zMin+zLift];
    
    % Their vertial concatenation is what you want
    pts = [P1; P2];
    
    % Because that's what line() wants to see
    line(pts(:,1), pts(:,2), pts(:,3),'Color',[1 0 1],'LineWidth',linWidth);
end

for i = 1:(numBlocks_M-1) % horizontal lines (along x)
    
    % Your two points
    P1 = [xyMin(1),yVec(i),zMin+zLift];
    P2 = [xyMax(1),yVec(i),zMin+zLift];
    
    % Their vertial concatenation is what you want
    pts = [P1; P2];
    
    % Because that's what line() wants to see
    line(pts(:,1), pts(:,2), pts(:,3),'Color',[1 0 1],'LineWidth',linWidth);
    
end

zLiftSource = 10;

scatter3(floodSource_xyz(:,1),floodSource_xyz(:,2),...
    floodSource_xyz(:,3)+zLiftSource,sourceMarkerSize,[1 0 0],'filled');
xlabel('x[m]');ylabel('y[m]');zlabel('z[m]');
title('Source points');

c = colorbar;
c.Label.String = 'z[m]';
caxis([cmin1,cmax1]);
view(0,90);

return
