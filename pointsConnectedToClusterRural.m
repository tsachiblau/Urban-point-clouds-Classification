function [newTmpLables] = pointsConnectedToCluster(xyz,pointLables,...
                          numOfNighbors,numOfTopNighbors,curvature)
    newTmpLables=pointLables;
    lengthOfPointLables=length(pointLables);
    tagNum=max(pointLables);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %print before adding
%     close all;
%     cloudColor= colorPointCloud(pointCloud(xyz),pointLables);
%     figure;
%     pcshow(cloudColor,'MarkerSize' ,40);
%     xlabel('x');ylabel('y');zlabel('z');
%     title('before graph');
%     daspect([1 1 1]);
% %   	return;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% build tree
    
    % if the group is empty
    if size(xyz,1)==0 | size(xyz,2)==0
        return;
    end
    [nighbors,distanceFromNighbors] = knnsearch(xyz,xyz,'k',numOfNighbors+1);
    
    %remove self distance
    nighbors(:,1)=[];
    distanceFromNighbors(:,1)=[];
    
    %get close clusteres
    closeDistanceFromNighborsIdx=distanceFromNighbors>0;
    s=[];
    t=[];
    w=[];

    for i=1:lengthOfPointLables
  
         tmp_t=nighbors(i,closeDistanceFromNighborsIdx(i,:))';
%         if i ==18490
%             x=1;
%         end
        
        %remove double edges
        if length(tmp_t)>0
            for j=1:length(tmp_t)
                if tmp_t(j)>i
                    %find the idx of the node in the new node
                    idxToRemove=find( nighbors(tmp_t(j),:)==i);
                    %change to not neighbors
                    closeDistanceFromNighborsIdx(tmp_t(j),idxToRemove)=false;
                end 
            end 
        end
        
        %check
%         if closeDistanceFromNighborsIdx(i,:)==1913
%             x=1;
%         end
        t=[t;tmp_t];
        s=[s;ones(length(tmp_t),1)*i];
        w=[w;distanceFromNighbors(i,closeDistanceFromNighborsIdx(i,:) )'];
    end
    
%%
    %build tree
    
    G=graph(s,t,w);
    %T=minspantree(G,'Type','forest');
    T=G;
    if size(T.Edges,1)>0
        maxDistance=max(T.Edges.Weight);
        f=@(x) 1./(x);
        T.Edges.Weight=f(T.Edges.Weight);
    end 
    
%     figure;
%     p1=plot(T,'XData',xyz(1:size(T.Nodes,1),1),'YData',xyz(1:size(T.Nodes,1),2),'ZData',xyz(1:size(T.Nodes,1),3));
% 
 %% connect unconected components
    bins = conncomp(T);
    uniqueBins=unique(bins);
    
    
    
    numOfMaxGroup=0;
    maxGroupQuantity=-1;
    % get max group
    for i=uniqueBins
        tmpQuantity=size(find(bins==i),2);
        if tmpQuantity>maxGroupQuantity
            maxGroupQuantity=tmpQuantity;
            numOfMaxGroup=i;
        end
%         close all;
%         figure;
%         p1=plot(T,'XData',xyz(1:size(T.Nodes,1),1),'YData',xyz(1:size(T.Nodes,1),2),'ZData',xyz(1:size(T.Nodes,1),3));
%         highlight(p1,find(bins==i),'NodeColor','green','MarkerSize',8);
    end
%     return;
    while size(uniqueBins,2)>1 & size(uniqueBins,1)>=1
        %x=size(uniqueBins,2)
        numOfseccondBigGroup=uniqueBins(find(uniqueBins~=numOfMaxGroup,1));
        bigGroupIdx=find(bins==numOfMaxGroup);
        seccondBigGroupIdx=find(bins==numOfseccondBigGroup);
        
        %add one weight between two groups
        [pointsToConnect,distanceOfpointsToConnect] = knnsearch(...
                                                                xyz(bigGroupIdx,:),...
                                                                xyz(seccondBigGroupIdx,:),'k',1);
        
        %get the idx of the closetst point                                                    
        if size(seccondBigGroupIdx,2)>1 & size(seccondBigGroupIdx,1)>=1
            [~,I]=min(distanceOfpointsToConnect,[],1);
            mainClusterClosestIdx=bigGroupIdx(pointsToConnect(I));
            otherClusterClosestIdx=seccondBigGroupIdx(1);
            tmpWeight=distanceOfpointsToConnect(I);
        else
            mainClusterClosestIdx=bigGroupIdx(pointsToConnect(1));
            otherClusterClosestIdx=seccondBigGroupIdx(1);
            tmpWeight=distanceOfpointsToConnect(1);
        end
        
        %pointsToConnect is the idx of the original pointLables vector
        %we need to add edge to the graph                                           
                                         
        T = addedge(T,mainClusterClosestIdx,otherClusterClosestIdx,f(tmpWeight));
        
        
        
        %%%%%%%%%%%%%%%%        %prepare for next iteration   %%%%%%%%%%%%%%%5
        
        %connect unconected components
        bins = conncomp(T);
        uniqueBins=unique(bins);

        numOfMaxGroup=0;
        maxGroupQuantity=-1;
        % get max group
        for i=uniqueBins
            tmpQuantity=size(find(bins==i),2);
            if tmpQuantity>maxGroupQuantity
                maxGroupQuantity=tmpQuantity;
                numOfMaxGroup=i;
            end
        end
    end
    
%     figure;
%     p1=plot(T,'XData',xyz(1:size(T.Nodes,1),1),'YData',xyz(1:size(T.Nodes,1),2),'ZData',xyz(1:size(T.Nodes,1),3));
    

%% get ready
    
 
    %get the cluster idx
    originalSource=find(~isnan(newTmpLables));
    source=find(~isnan(newTmpLables));
%% calc flowTH
    %calc flowTH
    sumOfmf=0;
    numberOftimes=1000;
    lengthOfSource=size(source,1);
    times=min(numberOftimes,floor(lengthOfSource/4));
    seeds=floor(linspace(1,lengthOfSource,times*2));
    mfVec=zeros(times,1);
%     x=0;
    for i=1:2:times*2
%         x=x+1;
%         if x==48
%             x=x;
%             %highlight(p1,191,'NodeColor','green','MarkerSize',8);
%         end
      	[mfVec(ceil(i/2)),~,~,~]=maxflow(T,source(seeds(i)),source(seeds(i+1)));  
        
    end
    
%     figure;
%     hist(mfVec);
    %flowTH=mean(mfVec)
%     flowTH=mean(T.Edges.Weight)+std(T.Edges.Weight);
    flowTH=mean(mfVec)

%     return;
%% calc curvetureTH
    tmpCurvature=curvature(source);
    curvetureMean=mean(tmpCurvature);
    goodCurv=tmpCurvature<curvetureMean;
    curvetureTH=mean(tmpCurvature(goodCurv));
    
    
%% add nans
    %get the nan idx
    sink=find(isnan(newTmpLables));
    
    %calc the min z of cluster
    minZofCluster=min(xyz(source,3));
%     figure;
%     p1=plot(T,'XData',xyz(1:size(T.Nodes,1),1),'YData',xyz(1:size(T.Nodes,1),2),'ZData',xyz(1:size(T.Nodes,1),3));
    partOfCluster=false;
    
   
    
    %if there is no sources
    if size(originalSource,1)==0 | size(originalSource,2)==0
        return;
    end
    
    %get nan curveture parameter
%     x=0;
    while size(sink,1)>0 & size(sink,2)>0
%         x=x+1;
%         if mod(x,100)==0
%             y=size(sink,1)
%             x=x
%         end
%         tmptmp=find(sink==7327);
%         if size(tmptmp,1)==0
%             x=x;
%         end
%         if x==108
%             x=x;
%         end
        [mf,~,cs,ct]=maxflow(T,originalSource(1),sink(1));
        
        %get better curvature
        rmLabledIdx=ismember(ct,originalSource);
        ctIdx=find(rmLabledIdx);
        if rmLabledIdx==0
            tmpCurv=mean(curvature(ct));
        else
            tmpCurv=curvature(sink(1));
        end
        
        %update weights
        mf=mf+getWeightToPoint(xyz(sink(1),:),flowTH,minZofCluster,curvature(sink(1)),curvetureTH,mf);
%         close all;
%         figure;
%         p1=plot(T,'XData',xyz(1:size(T.Nodes,1),1),'YData',xyz(1:size(T.Nodes,1),2),'ZData',xyz(1:size(T.Nodes,1),3));
%         highlight(p1,cs,'NodeColor','green','MarkerSize',9);
%         highlight(p1,3764,'NodeColor','green','MarkerSize',8);

%         %if its not connected
%         if mf==0
%             sink(1)=[];
%             continue;
%         end
%                
%         
        if mf>flowTH
            %change lable
            newTmpLables(sink(1))=tagNum;
            
            %remove from sink
            sink(1)=[];
            
            %update source
            source=find(~isnan(newTmpLables));
        else  %there is weak conectivity
            
            %there is more than 1 in the group
            if length(ct)>=1 & length(cs)>=1
                weChecked=[];
                allNan=false;
                while mf<flowTH & allNan==false & length(cs)>=1 & length(ct)>=1
                    %is nan?
                    rmLabledIdx=ismember(ct,originalSource);
                    ctIdx=find(rmLabledIdx);
                    
                    %the tags that we want to check
                    sourcesIdxNeedToCheck=ct(ctIdx);
        
                    if rmLabledIdx==0
                        allNan=true;
                    else
                        %take only the closest sorces
                        
                        [topListSourcesNeedToCheck,~] = knnsearch(xyz(sourcesIdxNeedToCheck,:),xyz(sink(1),:),'k',numOfTopNighbors);
                        
                        
                        topListSourcesNeedToCheckIdx=sourcesIdxNeedToCheck(topListSourcesNeedToCheck);
                        %check if we are in loop
                        looping= ismember(topListSourcesNeedToCheckIdx,weChecked);
                        if looping==1
                            break;
                        end
                        
                        goodIdx=find(~ismember(topListSourcesNeedToCheckIdx,weChecked),1);
                        weChecked=[weChecked,topListSourcesNeedToCheckIdx(goodIdx)];
                        if length(goodIdx)>0
        [mf,~,cs,ct]=maxflow(T,originalSource(1),sink(1));
        
        %get better curvature
        rmLabledIdx=ismember(ct,originalSource);
        ctIdx=find(rmLabledIdx);
        if rmLabledIdx==0
            tmpCurv=mean(curvature(ct));
        else
            tmpCurv=curvature(sink(1));
        end
        
        %update weights
        mf=mf+getWeightToPoint(xyz(sink(1),:),flowTH,minZofCluster,curvature(sink(1)),curvetureTH,mf);
%                             close all;
%                             figure;
%                             p1=plot(T,'XData',xyz(1:size(T.Nodes,1),1),'YData',xyz(1:size(T.Nodes,1),2),'ZData',xyz(1:size(T.Nodes,1),3));
%                             highlight(p1,cs,'NodeColor','green','MarkerSize',9);
%                             highlight(p1,ct,'NodeColor','red','MarkerSize',9);
                        end
                    end      
                end
                                
                %why wee out of the loop?
                if mf>flowTH  % we need to add it
                    %change lable
                    newTmpLables(sink(1))=tagNum;

                    %remove from sink
                    sink(1)=[];
            
                    %update source
                    source=find(~isnan(newTmpLables));
                
                elseif allNan==true  %remove all nan in ct  
                    %remove from sink
                    idxToRm=ismember(sink,ct);
                    sink(idxToRm)=[];
                else
                    %remove from sink
                    sink(1)=[];
                end
                
            else %we need to disconnect the sink
                %remove from sink
                sink(1)=[];
            end
            
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %print after adding
%     cloudColor= colorPointCloud(pointCloud(xyz),newTmpLables);
%     figure;
%     pcshow(cloudColor,'MarkerSize' ,40);
%     xlabel('x');ylabel('y');zlabel('z');
%     title('end graph');
%     daspect([1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end