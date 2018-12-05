function [newLables] = pumpShapes(lables,numOfPixels)
    
   
    %make copy of lables
    newLables=lables;
    
    
    [D,IDX] = bwdist(lables,'euclidean');
    
    %get all the relevant pixels
    notThelightened=D>0;
    closeEnough=D<=numOfPixels;
    needTolighten=notThelightened&closeEnough;
    for row=1:size(D,1)
        for column=1:size(D,2)
            if needTolighten(row,column)==1
               newLables(row,column)= lables(IDX(row,column));
            end
        end
    end

end