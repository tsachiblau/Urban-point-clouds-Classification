function [weight] = getWeightToPoint(pointXYZ,flowTH,minZofCluster,curvature,curvetureTH,mf)
    if mf<flowTH
        hightTH=2;
        certainBuildingsHeight=6;
        Hweight=0;
        x=(pointXYZ(3)-minZofCluster);
        if (x-hightTH)-certainBuildingsHeight>0
            Hweight=5*flowTH;
        elseif (x-hightTH)>0
            Hweight=1*(flowTH-mf)*(1/(certainBuildingsHeight-hightTH))*(x-hightTH);
        end

        Cweight=0;
        curvetureCoef=4;
        if curvature<curvetureTH
            tmp=curvetureTH-curvature;
            x=tmp/curvetureTH;
        else
            tmp=curvature-curvetureTH;
            x=-tmp/curvetureTH;
        end
        x=2*(1/(1+exp(-curvetureCoef*x))-0.5);
        Cweight=0.9*(flowTH-mf)*x;

        %sum the weights
        weight=Cweight+Hweight;
    else
        weight=0;
    end
end