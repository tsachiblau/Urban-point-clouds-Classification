function [weight] = getWeightToPoint(pointXYZ,flowTH,minZofCluster,curvature,curvetureTH,mf)
    
    hightTH=2;
    certainBuildingsHeight=6;
    Hweight=0;
    x=(pointXYZ(3)-minZofCluster);
    if (x-hightTH)-certainBuildingsHeight>0
        Hweight=5;
    elseif (x-hightTH)>0
        Hweight=2*(flowTH-mf)*(1/(certainBuildingsHeight-hightTH))*(x-hightTH);
    end
    
    Cweight=0;
    curvetureCoef=1;
    if curvature<curvetureTH
        tmp=curvetureTH-curvature;
        x=tmp/curvetureTH;
    else
        tmp=curvature-curvetureTH;
        x=-tmp/curvetureTH;
    end
    x=1.5*(1/(1+exp(-curvetureCoef*x))-0.5);
    Cweight=2*(flowTH-mf)*x;

    weight=Cweight+Hweight;

        
end