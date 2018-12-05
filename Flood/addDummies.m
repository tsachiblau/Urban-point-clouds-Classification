function [ xyzWithDummies] = addDummies( xyz , jump)
xyzWithDummies = [];

for i=1:(size(xyz,1)-1)
    dist = pdist( [ xyz(i,:) ; xyz((i+1),:) ],'euclidean');
    if dist > jump
        numOfPartitions = floor(dist/jump)+1;
        xArr = linspace(xyz(i,1),xyz(i+1,1),numOfPartitions)';
        yArr = linspace(xyz(i,2),xyz(i+1,2),numOfPartitions)';
        zArr = linspace(xyz(i,3),xyz(i+1,3),numOfPartitions)';
        xyzWithDummies = [xyzWithDummies ; [xArr , yArr , zArr ]];
    else
        xyzWithDummies = [ xyzWithDummies ; [xyz(i,:);xyz(i+1,:)]];
    end
end

end

