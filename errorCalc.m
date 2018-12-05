function errorCalc( ourLables , realLables ,inputFlag, origCloud)

if inputFlag == 0
    groundIdx = realLables==202030000 | realLables==202000000 | realLables==202010000 | realLables==202020000 | realLables==202040000;
    realLables(groundIdx)=1;
    
    buildingIdx = realLables==203000000;
    realLables(buildingIdx)=2;
    
    untagsIdx= realLables~=1 & realLables~=2;
    realLables(untagsIdx)=0;
elseif inputFlag == 1
    realLables(realLables == 3) = 0;
elseif inputFlag == 2
    realLables(realLables==1 | realLables ==2) = 1;
    realLables(realLables==3 | realLables ==4 | realLables ==6| realLables ==7| realLables ==8) = 0;
    realLables(realLables ==5) = 2;
end

% real
realLables_OthersLogic= realLables==0;      % vector of '1' where there is "others" in the realTags vector
realLables_GroundLogic= realLables==1;      % vector of '1' where there is "ground" in the realTags vector
realLables_BuildingsLogic= realLables==2;   % vector of '1' where there is "building" in the realTags vector
% our
ourLables_OthersLogic= ourLables==0;       % vector of '1' where there is "others" in our Tags vector
ourLables_GroundLogic= ourLables==1;      % vector of '1' where there is "ground" in our Tags vector
ourLables_buildingsLogic= ourLables==2;    % vector of '1' where there is "building" our Tags vector

%% comapare ground detection
TP_logical_G = realLables_GroundLogic & ourLables_GroundLogic;
TN_logical_G = (~realLables_GroundLogic) & (~ourLables_GroundLogic);
FN_logical_G = realLables_GroundLogic & (~ourLables_GroundLogic);
FP_logical_G = (~realLables_GroundLogic) & ourLables_GroundLogic;

Ground_TP= sum( double( TP_logical_G));
Ground_TN= sum( double( TN_logical_G));
Ground_FN= sum( double( FN_logical_G));
Ground_FP= sum( double( FP_logical_G));

groundRec=( Ground_TP / (Ground_TP + Ground_FN) );
display(['ground recall messure is:  ',num2str(groundRec)]);

groundPre=( Ground_TP / (Ground_TP+Ground_FP) );
display(['ground Precision  messure is:  ',num2str(groundPre)]);

groundfMeasure=2*(groundRec*groundPre)/(groundRec+groundPre);
display(['ground F-messure is:  ',num2str(groundfMeasure)]);


%% comapare "Others" detection

TP_logical_O = realLables_OthersLogic & ourLables_OthersLogic;
TN_logical_O = (~realLables_OthersLogic) & (~ourLables_OthersLogic);
FN_logical_O = realLables_OthersLogic & (~ourLables_OthersLogic);
FP_logical_O = (~realLables_OthersLogic) & ourLables_OthersLogic;

Others_TP= sum( double( TP_logical_O));
Others_TN= sum( double( TN_logical_O));
Others_FN= sum( double( FN_logical_O));
Others_FP= sum( double( FP_logical_O));


OthersRec=( Others_TP / (Others_TP + Others_FN) );
display(['Others recall messure is:  ',num2str(OthersRec)]);

OthersPre=( Others_TP / (Others_TP + Others_FP) );
display(['Others Precision  messure is:  ',num2str(OthersPre)]);

OthersfMeasure=2*(OthersRec*OthersPre)/(OthersRec+OthersPre);
display(['Others F-messure is:  ',num2str(OthersfMeasure)]);

%% comapare "building" detection

TP_logical_B = realLables_BuildingsLogic & ourLables_buildingsLogic;
TN_logical_B = (~realLables_BuildingsLogic) & (~ourLables_buildingsLogic);
FN_logical_B = realLables_BuildingsLogic & (~ourLables_buildingsLogic);
FP_logical_B = (~realLables_BuildingsLogic) & ourLables_buildingsLogic;

Buildings_TP= sum( double( TP_logical_B));
Buildings_TN= sum( double( TN_logical_B));
Buildings_FN= sum( double( FN_logical_B));
Buildings_FP= sum( double( FP_logical_B));

BuildingsRec=( Buildings_TP / (Buildings_TP + Buildings_FN) );
display(['Buildings recall messure is:  ',num2str(BuildingsRec)]);

BuildingsPre=( Buildings_TP / (Buildings_TP+Buildings_FP) );
display(['Buildings Precision  messure is:  ',num2str(BuildingsPre)]);

BuildingsfMeasure=2*(BuildingsRec*BuildingsPre)/(BuildingsRec+BuildingsPre);
display(['Buildings F-messure is:  ',num2str(BuildingsfMeasure)]);

% %% show
    % ground
% 
%     combinedLables=double(TP_logical_G);
%     combinedLables(TP_logical_G)=0;
%     combinedLables(TN_logical_G)=1;
%     combinedLables(FN_logical_G)=2;
%     combinedLables(FP_logical_G)=3;
% % green = TP, red = FN,TN, blue = FP
%     figure;
%     cloud_colored = colorPointCloud(origCloud,combinedLables);
%     pcshow(cloud_colored);
% 
%     % buildings
%     combinedLables=double(TP_logical_B);
%     combinedLables(TP_logical_B)=0;
%     combinedLables(TN_logical_B)=1;
%     combinedLables(FN_logical_B)=2;
%     combinedLables(FP_logical_B)=3;
%     figure;
%     cloud_colored = displayPointCloudLabels(DS_cloud,combinedLables);
%     pcshow(cloud_colored);
confusionMatrix = zeros(3);
confusionMatrix(1,1) = sum(ourLables_OthersLogic & realLables_OthersLogic);
confusionMatrix(1,2) = sum(ourLables_GroundLogic & realLables_OthersLogic);
confusionMatrix(1,3) = sum(ourLables_buildingsLogic & realLables_OthersLogic);
confusionMatrix(2,1) = sum(ourLables_OthersLogic & realLables_GroundLogic);
confusionMatrix(2,2) = sum(ourLables_GroundLogic & realLables_GroundLogic);
confusionMatrix(2,3) = sum(ourLables_buildingsLogic & realLables_GroundLogic);
confusionMatrix(3,1) = sum(ourLables_OthersLogic & realLables_BuildingsLogic);
confusionMatrix(3,2) = sum(ourLables_GroundLogic & realLables_BuildingsLogic);
confusionMatrix(3,3) = sum(ourLables_buildingsLogic & realLables_BuildingsLogic);

    overall_accuracy = sum(double(eq(realLables,ourLables))) / length(realLables);
    display(['overall_accuracy is:  ',num2str(overall_accuracy)]);

end



