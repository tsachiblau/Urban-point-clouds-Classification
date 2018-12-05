function [minMap,maxMap] = ...
    fillPixels_NN(minMap,maxMap,dataMask,pixels2fill)

%% Find closest non-zero (non-empty) pixel for each pixel
[~,IDX] = bwdist(dataMask);

%% Fill in empty pixels
minMap(pixels2fill) = minMap(IDX(pixels2fill));
maxMap(pixels2fill) = maxMap(IDX(pixels2fill));
% imshow(minMap,[])
% imshow(maxMap,[])

return
