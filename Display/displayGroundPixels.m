function displayGroundPixels(groundBinary,groundMask,heightMap,cloud)

defaultFontSize = 20;
set(0,'DefaultAxesFontSize',defaultFontSize);
set(0, 'DefaultAxesFontName','Times');

xyMin = [cloud.XLimits(1) cloud.YLimits(1)];
xyMax = [cloud.XLimits(2) cloud.YLimits(2)];
xLimDisp = [xyMax(1); xyMin(1)];
yLimDisp = [xyMax(2); xyMin(2)];

%%
figure;
imagesc(xLimDisp,yLimDisp,rot90(groundBinary,2));
set(gca,'YDir','normal');
xlabel('x[m]');ylabel('y[m]');
title(['Ground mask']);
colormap gray


%%
groundBinary_disp = uint8(255*(groundMask<inf));
groundBinary_disp(heightMap==inf) = 128;

figure;
imagesc(xLimDisp,yLimDisp,rot90(groundBinary_disp,2));
set(gca,'YDir','normal');
xlabel('x[m]');ylabel('y[m]');
colormap gray
colorbar;
colorbar('Ticks',[0,128,255],'TickLabels',{'Not Ground','No Data','Ground'});
title('Ground Segmentation Result');

%% height colored segmentation
heightMap_seg = heightMap;
heightMap_seg(not(groundBinary)) = inf;

figure;
imagesc(xLimDisp,yLimDisp,rot90(heightMap_seg,2));
set(gca,'YDir','normal');
xlabel('x[m]');ylabel('y[m]');
colormap jet
c = colorbar;
[cmin1 cmax1] = caxis;
c.Label.String = 'z[m]';
title('Height-map flood result');

%%
heightMap_seg2 = heightMap;
heightMap_seg2(groundBinary) = inf;

figure;
imagesc(xLimDisp,yLimDisp,rot90(heightMap_seg2,2));
set(gca,'YDir','normal');
xlabel('x[m]');ylabel('y[m]');
colormap jet
c = colorbar;
% caxis([cmin1 cmax1]);
% caxis([0 80]);
c.Label.String = '\fontsize{16}z[m]';
title('"Not ground" Segmentation Result - height colored');



return
