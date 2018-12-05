function [floodSource,numBlocks_M,numBlocks_N] = ...
    findFloodSources(minMap,differenceMap,...
    minBlockSide,maxZDiff)

%% Find source pixels
possibleSourceMap = minMap;
possibleSourceMap(differenceMap > maxZDiff) = nan;
% imshow(possibleSourceMap,[])

% determine number of blocks
[M,N] = size(minMap);

numBlocks_M = floor(M/minBlockSide);
numBlocks_N = floor(N/minBlockSide);

% divide cloud to blocks
blocks_cell = cell(numBlocks_M,numBlocks_N);
blockHeight = floor(M/numBlocks_M);
blockWidth = floor(N/numBlocks_N);

for i = 1:numBlocks_M
    m_Range = [1:blockHeight] + (i-1)*blockHeight;
    for j = 1:numBlocks_N
        n_Range = [1:blockWidth] + (j-1)*blockWidth;
        blocks_cell{i,j}.data = possibleSourceMap(m_Range,n_Range);
        blocks_cell{i,j}.topLeftCorner = [(1 + (i-1)*blockHeight) (1 + (j-1)*blockWidth)];
        
        minMapBlock = minMap(m_Range,n_Range);
        % imshow(minMapBlock,[])
        blocks_cell{i,j}.numValidPixels = sum(minMapBlock(:)<inf);        
    end
end

% find min pixel in each block
findMinPixelInBlock_handle = @(x) findMinPixelInBlock(x,minBlockSide,1);

[minPixelInds,blockType,~] = ...
    cellfun(findMinPixelInBlock_handle,blocks_cell,'UniformOutput',false);

blockType = cell2mat(blockType);
blockValid = (blockType==2);

%%
minPixelInds = minPixelInds(:);
floodSource = cell2mat(minPixelInds);

floodSource_isValid = blockValid(:);
floodSource(not(floodSource_isValid),:) = [];

return
