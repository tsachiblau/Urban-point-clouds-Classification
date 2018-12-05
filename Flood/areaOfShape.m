function [ area ] = areaOfShape( xyz,R )
    projection=xyz*R;
    max_of_profjection=max(projection);
    min_of_profjection=min(projection);
    area=(max_of_profjection(1)-min_of_profjection(1))*(max_of_profjection(2)-min_of_profjection(2));
end

