function displayHeightMaps(minMap,maxMap,differenceMap,cloud)

%% Display
xyMin = [cloud.XLimits(1) cloud.YLimits(1)];
xyMax = [cloud.XLimits(2) cloud.YLimits(2)];
xLimDisp = [xyMax(1); xyMin(1)];
yLimDisp = [xyMax(2); xyMin(2)];

[M,N] = size(minMap);

defaultFontSize = 16;
set(0,'DefaultAxesFontSize',defaultFontSize);
set(0, 'DefaultAxesFontName','Times');

%% Display min-map
figure;
imagesc(xLimDisp,yLimDisp,rot90(minMap,2));
set(gca,'YDir','normal');
colormap jet
xlabel('x[m]');ylabel('y[m]');
title(['Min. height map - ',num2str(M),'x',num2str(N)]);
c = colorbar;
c.Label.String = 'z[m]';
[cmin1,cmax1] = caxis;
axis equal

%% Display max-map
figure;
imagesc(xLimDisp,yLimDisp,rot90(maxMap,2));
set(gca,'YDir','normal');
colormap jet
xlabel('x[m]');ylabel('y[m]');
title(['Max. height map - ',num2str(M),'x',num2str(N)]);
c = colorbar;
c.Label.String = 'z[m]';
caxis([cmin1 cmax1]);
axis equal

%% Display max-min
figure;
imagesc(xLimDisp,yLimDisp,rot90(differenceMap,2));
set(gca,'YDir','normal');
colormap jet
xlabel('x[m]');ylabel('y[m]');
title(['Max-Min height map - ',num2str(M),'x',num2str(N)]);
c = colorbar;
c.Label.String = 'Max(z) - Min(z) [m]';
axis equal

return
