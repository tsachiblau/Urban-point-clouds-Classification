function displayCurv(xyz,curveture)
    
    %cloud = pointCloud(xyz);
    
    

    %print curvature
    
    figure;
    pcshow(xyz,curveture);
    xlabel('x');ylabel('y');zlabel('z');
    c = colorbar;
    c.Label.String = 'Curvature';
    title('Estimated curvature');
    colormap jet
    caxis([0 1/3]);
    daspect([1 1 1]);


end