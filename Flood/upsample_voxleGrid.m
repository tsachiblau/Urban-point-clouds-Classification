function [ finalTags ] = upsample_voxleGrid( groundCompleteCloudIdx, nonGroundCompleteCloudIdx ,OurDSTags,ind )

finalTags = nan(length(groundCompleteCloudIdx),1);
finalTags(groundCompleteCloudIdx) = 1;
tmpTags = zeros(sum(single(nonGroundCompleteCloudIdx)),1);
tmpTags = OurDSTags(ind);
finalTags(nonGroundCompleteCloudIdx) = tmpTags;

end

