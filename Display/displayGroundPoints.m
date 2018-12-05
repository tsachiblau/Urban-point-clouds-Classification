function displayGroundPoints(cloud,groundPointsFlags)

xyz = cloud.Location;
xyzGround = xyz(groundPointsFlags,:);
xyzNotGround = xyz(not(groundPointsFlags),:);

%% Display non-ground points
figure;
pcshow(xyzNotGround);
xlabel('x[m]');ylabel('y[m]');zlabel('z[m]');
c = colorbar;
c.Label.String = 'z[m]';
colormap jet
title('Non-ground points');
view(0,90);

xlim1 = xlim;ylim1 = ylim;

%% Display ground points
figure;
pcshow(xyzGround);
xlabel('x[m]');ylabel('y[m]');zlabel('z[m]');
c = colorbar;
c.Label.String = 'z[m]';
colormap jet
title('Ground points');
view(0,90);

xlim(xlim1);ylim(ylim1);

return
