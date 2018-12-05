function [minPixelInds,blockType,minPixelHeight] = ...
    findMinPixelInBlock(block,minBlockSide,C)

data = block.data;
[M,N] = size(data);

% imshow(data,[])

[minPixelHeight,linInd] = min(data(:));
[m,n] = ind2sub([M,N],linInd);
minPixelInds = [m-1 n-1] + block.topLeftCorner;

numValidPixels = block.numValidPixels;

if numValidPixels == 0
    blockType = 0; % empty

elseif numValidPixels < C*(minBlockSide^2)
    blockType = 1; % partial
else
    blockType = 2; % valid
end

return
