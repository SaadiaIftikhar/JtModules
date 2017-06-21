function CutMask=PerimeterSegmentation(LabelImage,PerimeterTrace,MaxEqivRadius,MinEquivAngle,AngleDeviation,MaxDistance,DistanceAngleRatio,DistanceMethod,AngleMethod,CutMethod,MinArea)
%PERIMETERSEGMENTATION segment objects based on shape
%   CUTMASK=PERIMETERSEGMENTATION(LABELIMAGE,PERIMETERTRACE,MAXEQUIVRADIUS,MINEQUIVANGLE,ANGLEDEVIATION,MAXDISTANCE,DISTANCEANGLERATIO,CUTPOINT)
%   seperates clumped objects (native shape round) based on a dumbell like
%   appereance. Clumped objects are characterized by opposite concave
%   regions.
%
%   CUTMASK returns a logical image the size of LABELIMAGE indicating the
%   cuts necessary to segment the objects (cut=1 background=0)
%
%   LABELIMAGE as created by bwlabel represents the objects
%
%   PERIMETERTRACE is a cell array as created by PerimeterAnalysis.m
%   containing the properties of the object perimeters.
%
%   PERIMETERSEGMENTATION cuts objects based on concave regions of objects:
%   Concave region below a certain equivalent radius (MAXEQUIVRADIUS in px)
%   and above an equivalent circle segment (MINEQUIVANGLE in radian)
%   qualify for cutting.
%
%   All possible pairs of qualifying regions are scored by distance and
%   direction of the curvature. Ideally two opposing concave regions have
%   a 180 degree angle between their curvature vectors.
%
%   ANGLEDEVIATION gives the maximum deviation (in radians) from the ideal
%   geometry. The angle measure depends on the setting of
%   ANGLEMETHOD (angle scoring):
%       'center' angle between the mean curvature of the regions
%       'curvature' angle between the normal vectors at the points of 
%       	maximum curvature of the regions
%       'best' best possible angle between a pair of normal vectors of
%       	the regions
%       'best_inline' pick pair of normal vectors minimizing the mean angle
%       	to the cut line between regions. ANGLEDEVIATION values should
%           be halfed for this method (90 degree deviation for each vector
%           from the connecting line makes the vectors parallel !)  
%       	                  
%   MAXDISTANCE gives the maximum distance between the regions. The
%   distance measure depends on the setting of
%   DISTANCEMETHOD (distance scoring):
%       'center' distance metween the centers of the regions
%       'curvature' distance between the maximum curvature points
%                   of the regions
%       'best' minimum distance between the regions
%
%   DISTANCEANGLERATIO gives the weighting ratio for calculating the final
%   scoring matrix (1=pure angle, 0=pure distance score). Based on this
%   score mutually closest pairs are selected.
%
%   Currently only pairwise cuts are performed (only one cut per region).
%   In some cases it may be beneficial to apply the function two times in
%   a row in order to perform multiple cuts in one concave region.
%
%   CUTMETHOD defines the way the cut is placed between the regions:
%       'angle' cuts according to the ANGLEMETHOD
%       'distance' cuts according to the DISTANCEMETHOD
%       Influence of DISTANCEMETHOD/ANGLEMETHOD:
%           'angle' or 'distance'+'center' cuts between the centers of the
%               two regions
%           'angle' or 'distance'+'curvature' cuts between the maximum 
%               curvature points of the
%                       regions
%           'distance'+'best' cuts at the minimum distance/best angle 
%               alignment between the regions
%           'angle'+'best' cuts between the points of the best aligning
%               normal vectors
%           'angle'+'best_inline' cuts between the points where the normal
%               vectors have the smallest angle to the cutting line (normal
%               vectors as parralel as possible to the cutting line)                      
%
%   MINAREA prevents single cuts from producing too small objects. A cut is
%   prohibited if one of the resulting parts is smaller than MINAREA.
%   Nevertheless multiple cuts can still produce objects smaller than
%   MINAREA.
%
%   [Anatol Schwab 21.11.2012]

%load('temp_img.mat');
%LoGLabel=bwlabel(LoGMask);
%testv=PerimeterAnalysis(LoGLabel,10);
%testw=PerimeterSegmentation(LoGLabel,testv,50,0.1*2*pi,0.2*2*pi,50,0.5,'best','best','angle',500);
%testwB=PerimeterSegmentation(LoGLabelB,testvB,50,0.05*2*pi,0.2*2*pi,50,0.5,'distance',500);%second pass
%TODO: AngleScore scaling  should depend on AngleDeviation, Distance Score
%scaling.....?
%Reapplying sometimes helps, because no tri-cuts (one region has two cuts
%to two other regions) are considered!
%TODO: variable cut line width, fix issue with connecting pixels at end
%TODO: maybe allow mono-cuts for regions with extreme high curvature?

ImSize=size(LabelImage);
ObjectIDs=setdiff(unique(LabelImage(:)),0);%just to make sure
NumObjects=length(ObjectIDs);
CutCoordList=[];%contains lines([r1,c1,r2,c2;...]) to perform cuts, %Attention: CutList is recalculated for each object, CutCoordList is accumulated!

for i=1:NumObjects
    %% Load concave regions for current object
    CurrentObject=ObjectIDs(i);
    CurrentPreimProps=PerimeterTrace{i};
    ConcaveRegions=bwlabel((CurrentPreimProps(:,11)==-1));%bwlabel also works nicely in 1D!
    NumConcave=length(setdiff(unique(ConcaveRegions),0));
    ConcaveRegionProps=zeros(NumConcave,12,'double');%summary of properties for each concave region:
    %% Characterize concave regions of current object
    for j=1:NumConcave
        CurrentRegionProps=CurrentPreimProps(ConcaveRegions==j,:);
        Curvature=CurrentRegionProps(:,5:6);
        NormalVectors=CurrentRegionProps(:,3:4);
        NormCurvature=CurrentRegionProps(:,9);
        ConcaveRegionProps(j,1)=max(NormCurvature);
        ConcaveRegionProps(j,2)=mean(NormCurvature);
        MaximaIndices=(NormCurvature==ConcaveRegionProps(j,1));%indices of pixels where the curvature is maximal
        ConcaveRegionProps(j,3:4)=mean(NormalVectors(MaximaIndices,:),1);%normal vector at pixels where the norm of the curvature is maximal (the maximum may apper multiple times in few cases ->mean)
        ConcaveRegionProps(j,5:6)=mean(NormalVectors,1);%mean normal vector
        FirstMaximumIndex=find(MaximaIndices,1,'first');
        LastMaximumIndex=find(MaximaIndices,1,'last');
        MeanMaximumIndex=round((LastMaximumIndex+FirstMaximumIndex)/2);%index of the pixel at which the curvature is maximal, if there are more than one maxima, a pixel inbetween the maxima is picked.
        ConcaveRegionProps(j,7:8)=CurrentRegionProps(MeanMaximumIndex,1:2);%pixel coordinates of maximum curvature, or mean between maxima
        ConcaveRegionProps(j,9:10)=CurrentRegionProps(round((size(CurrentRegionProps,1)+1)/2),1:2);%center pixel
        %The sum of curvature has a maximum of 2*pi (object shape=circle)
        %exceptions are very small object sizes or very large objects in
        %combination with a too small sliding window (eg window<r/15)
        %the total curvature is the equvalent angle(radian) of the curve!!!
        ConcaveRegionProps(j,11)=sum(NormCurvature);%total curvature of region = EQUIVALENT ANGLE of circle segemnt
        ConcaveRegionProps(j,12)=length(NormCurvature)/sum(NormCurvature);%EQUIVALENT RADIUS
        ConcaveRegionProps(j,13)=j;%save region ID so we can later quickly refer to information in CurrentPreimProps (see start of loop)
    end
    %% Select regions meeting the Radius/Angle criteria, score them
    
    %filter based on equivalent radius/angle
    QualifyingRegionsMask=(ConcaveRegionProps(:,11)>=MinEquivAngle)&(ConcaveRegionProps(:,12)<=MaxEqivRadius);
    SelectedRegions=ConcaveRegionProps(QualifyingRegionsMask,:);
    NumSelected=size(SelectedRegions,1);
    
    %calculate angle between curvature of all selected regions
    AngleCutPoints=cell(NumSelected);%AngleCutPoints and DistanceCutPoints are the same for center and curvature
    AngleMatrix=zeros(NumSelected,'double');
    for j=1:NumSelected
        for k=1:j%loop over half of matrix
            if strcmp(AngleMethod,'center')%angle between mean gradient of regions
                CurrentAngle=acos(dot(SelectedRegions(j,5:6),SelectedRegions(k,5:6))/(norm(SelectedRegions(j,5:6))*norm(SelectedRegions(k,5:6))));%mean
                CurrentCut=[SelectedRegions(j,[9,10]),SelectedRegions(k,[9,10])];
            elseif strcmp(AngleMethod,'curvature')%angle between maximum gradient of regions
                CurrentAngle=acos(dot(SelectedRegions(j,3:4),SelectedRegions(k,3:4))/(norm(SelectedRegions(j,3:4))*norm(SelectedRegions(k,3:4))));%max   
                CurrentCut=[SelectedRegions(j,[7,8]),SelectedRegions(k,[7,8])];
            elseif strcmp(AngleMethod,'best')||strcmp(AngleMethod,'best_inline')%best angle between regions
                RegionA=CurrentPreimProps(ConcaveRegions==SelectedRegions(j,13),[1 2 3 4]);%load both regions: coordinates AND normal vectors!
                RegionB=CurrentPreimProps(ConcaveRegions==SelectedRegions(k,13),[1 2 3 4]);%load both regions: coordinates AND normal vectors!
                AllAngles=zeros(size(RegionA,1),size(RegionB,1),'double');%A=row, B=colums (opposite for distance calc!)
                for l=1:size(RegionA,1)%full loop over all combinations of cutting between region A and B!
                    for m=1:size(RegionB,1)
                        if strcmp(AngleMethod,'best')
                            %'best' only optimizes for 180 degree vector alignment but not that vectors actually point at each other!
                            AllAngles(l,m)=acos(dot(RegionA(l,3:4),RegionB(m,3:4))/(norm(RegionA(l,3:4))*norm(RegionB(m,3:4))));
                        elseif strcmp(AngleMethod,'best_inline')
                            %'best inline' optimizes for the smallest angle of two vectors to the connecting line!
                            ConnectingVectorAB=RegionB(m,1:2)-RegionA(l,1:2);%vector from pixel in region A to B
                            ConnectingVectorAB=ConnectingVectorAB/norm(ConnectingVectorAB);
                            ConnectingVectorBA=-ConnectingVectorAB;
                            AngleDeviationA=acos(dot(RegionA(l,3:4),ConnectingVectorAB));%Angle between the normal vector and the connectiong vector=deviation from ideal geometry
                            AngleDeviationB=acos(dot(RegionB(m,3:4),ConnectingVectorBA));%Angle between the normal vector and the connectiong vector=deviation from ideal geometry
                            MeanAngleDeviation=(AngleDeviationA+AngleDeviationB)/2;%ranges between 0 and 180 degree, 0 being the ideal geometry
                            AllAngles(l,m)=pi-MeanAngleDeviation;%this is done just to conform with the further angle scoring where 180 is considered best, 0 worst
                        end
                    end
                end
                AllAngles=real(AllAngles);%maybe not necessary
                CurrentAngle=max(AllAngles(:));%max(acos)=pi, exactly what we are looking for (180 degree geometry)
                [PosA,PosB]=find(AllAngles==CurrentAngle);
                CurrentCut=[RegionA(PosA(1),1:2),RegionB(PosB(1),1:2)];%take first occurence of maximum if multiple maxima
            else
                error('%s: Angle score method unspecified, use either ''center'', ''curvature'', ''best'' or ''best_inline''\n',mfilename);
            end         
            AngleMatrix(j,k)=real(CurrentAngle);
            AngleMatrix(k,j)=real(CurrentAngle);
            AngleCutPoints{j,k}=CurrentCut;
            AngleCutPoints{k,j}=[CurrentCut(3:4),CurrentCut(1:2)];%directional cut...
        end
    end
    AngleScoreMask=(AngleMatrix>=(pi-AngleDeviation));%blanks all antires that deviate too much from the ideal (180 degree=pi) geometry
    AngleScore=AngleMatrix.*AngleScoreMask;
    %AngleScore=AngleScore/pi;%values now range between 0 and 1 (1=best)
    AngleScore=(AngleScore-(pi-AngleDeviation)).*AngleScoreMask./AngleDeviation;%better score: scale ranges form 1 (best, 180 degree=pi) to 0 (worst, everything >=pi-AngleDeviation) 
    %diagonal is zero already because angle between self and self is 0!
    
    %calculate distances between all selected regions    
    DistanceMatrix=zeros(NumSelected,'double');%AngleCutPoints and DistanceCutPoints are the same for center and curvature
    DistanceCutPoints=cell(NumSelected);%already save all possible cut points 
    for j=1:NumSelected
        for k=1:j
            %CurrentDistance=norm(SelectedRegions(j,9:10)-SelectedRegions(k,9:10));%default?
            if strcmp(DistanceMethod,'center')%distance between centers of regions
                CurrentDistance=norm(SelectedRegions(j,9:10)-SelectedRegions(k,9:10));
                CurrentCut=[SelectedRegions(j,[9,10]),SelectedRegions(k,[9,10])];
            elseif strcmp(DistanceMethod,'curvature')%distance between points of maximum gradient
                CurrentDistance=norm(SelectedRegions(j,7:8)-SelectedRegions(k,7:8));
                CurrentCut=[SelectedRegions(j,[7,8]),SelectedRegions(k,[7,8])];
            elseif strcmp(DistanceMethod,'best')%minimum distance between regions
                RegionA=CurrentPreimProps(ConcaveRegions==SelectedRegions(j,13),1:2);%load both regions: coordinates ONLY
                RegionB=CurrentPreimProps(ConcaveRegions==SelectedRegions(k,13),1:2);%load both regions: coordinates ONLY
                [MeshARow,MeshBRow]=meshgrid(RegionA(:,1),RegionB(:,1));
                [MeshACol,MeshBCol]=meshgrid(RegionA(:,2),RegionB(:,2));
                RowDist=MeshARow-MeshBRow;
                ColDist=MeshACol-MeshBCol;
                %Region dist: RegionA: columns, RegionB: rows
                RegionDist=sqrt(RowDist.*RowDist+ColDist.*ColDist);%euclidian distance for all possible cuts.
                CurrentDistance=min(RegionDist(:));
                [PosB,PosA]=find(RegionDist==CurrentDistance);%Notice: PosB=row, PosA=column in RegionDist
                CurrentCut=[RegionA(PosA(1),1:2),RegionB(PosB(1),1:2)];
            else
                error('%s: Distance score method unspecified, use either ''center'', ''curvature'' or ''best''\n',mfilename);
            end
            %matrix currently symmetric, no directionality of cuts
            DistanceMatrix(j,k)=CurrentDistance;
            DistanceMatrix(k,j)=CurrentDistance;
            DistanceCutPoints{j,k}=CurrentCut;
            DistanceCutPoints{k,j}=[CurrentCut(3:4),CurrentCut(1:2)];%directional cut...
        end
    end
    DistanceScoreMask=(DistanceMatrix<MaxDistance).*~diag(diag(ones(size(DistanceMatrix,1))));%only allow cuts to a certain length, delete diagonal
    DistanceScore=DistanceMatrix.*DistanceScoreMask;
    if ~isempty(max(DistanceScore(:)))
        %create linear score: 1 = best (0px distance), 0=worst (distance>=MaxDistance)
        DistanceScore=1-(DistanceScore./MaxDistance);
        DistanceScore=DistanceScore.*DistanceScoreMask;
    end
    
    %calculate partitioning score (avoid small bits cut off, maybe favour
    %equal size cuts)
    PartitionScore=ones(NumSelected,'double');
    PartitionMatrix=ones([NumSelected,NumSelected,2],'double');%(:,:,1) always contains the area of the larger part, 'ones' to have diagonal nonzero and avoid division by zero when calculating relative area
    if strcmp(CutMethod,'distance')
        CutPoints=DistanceCutPoints;
    elseif strcmp(CutMethod,'angle')
        CutPoints=AngleCutPoints;
    end
    %Use the precalculated cut points to determine the partitioning.
    %Currently each cut is considered on its own. Therefore it does not
    %make sense to enforce equal partitioning as the object could be cut
    %into multiple pieces in the end. While cutting off small pieces is
    %effectively avoided, areas arising from two cuts may still be smaller
    %than MinArea
    for j=1:NumSelected
        for k=1:j-1%exclude diagonal.
            %find indices in object perimeter list
            CutPointIndices=sort(find(ismember(CurrentPreimProps(:,1:2),[CutPoints{j,k}(1:2);CutPoints{j,k}(3:4)],'rows')));%locate the indices of the cut points in the list. the result is unambiguous since each border pixel(row,col) may only occur once in the list!
            CutPointIndexA=CutPointIndices(1);%lower index
            CutPointIndexB=CutPointIndices(2);%higher index
            RegionA=CurrentPreimProps(CutPointIndexA:CutPointIndexB,1:2);%Region from A to B
            RegionB=[CurrentPreimProps(CutPointIndexB:end,1:2);CurrentPreimProps(1:CutPointIndexA,1:2)];%Region from B to A
            AreaA=polyarea(RegionA(:,1),RegionA(:,2));
            AreaB=polyarea(RegionB(:,1),RegionB(:,2));
            PartitionMatrix(j,k,1)=max([AreaA,AreaB]);
            PartitionMatrix(k,j,1)=PartitionMatrix(j,k,1);
            PartitionMatrix(j,k,2)=min([AreaA,AreaB]);
            PartitionMatrix(k,j,2)=PartitionMatrix(j,k,2);
            if PartitionMatrix(j,k,1)<MinArea || PartitionMatrix(j,k,2)<MinArea 
                PartitionScore(j,k)=0;
                PartitionScore(k,j)=0;
            end
        end
    end
    
    %calculate final score
    OverallScore=AngleScore*DistanceAngleRatio+DistanceScore*(1-DistanceAngleRatio);%allow to set weight of different scores
    OverallScore=OverallScore.*DistanceScoreMask.*AngleScoreMask;%ensure that if one of the criteria is not met at all (=0) the other can't compensate!
    OverallScore=OverallScore.*PartitionScore;%avoid cutting off small pieces!
    %% Select PAIRS for cutting based on score
    [SortedScore,SortOrder]=sort(OverallScore,2,'descend');%SortOrder(RegionIdx,n) gives the n-th nearest neighbour index for the region with RegionIdx. eg SortOrder(3,1) returns the index of the closest(1) region to the region with index 3
    %check if regions are mutual nearest, add them to the cut list
    CutList=[];
    %Attention: CutList is recalculated for each object, CutCoordList is accumulated!
    %Attention: CutList contains indices in SelectedRegions, NOT ConcaveRegionProps, NOT CurrentPreimProps!
    %To get index in ConcaveRegionProps use: SelectedRegions(idx,13)
    %To get indices in CurrentPreimProps use: CurrentPreimProps(ConcaveRegions==SelectedRegions(idx,13),:);
    for j=1:size(SortOrder,1)
        if SortOrder(SortOrder(j,1),1)==j && j<SortOrder(j,1) && SortedScore(j,1)>0 && SortedScore(SortOrder(j,1),1)>0%Only connect regions that mutually the clostest ones! Only add pair ONCE to list! (via < ), also check for nonzero score!
            CutList=[CutList;[j,SortOrder(j,1)]];%cut between j and SortOrder(j,1)
        end
    end
    
    %% Process CutList, generate cut endpoint coordinates
    for j=1:size(CutList,1)
        if strcmp(CutMethod,'distance')
            CutCoordList=[CutCoordList;DistanceCutPoints{CutList(j,1),CutList(j,2)}];
        elseif strcmp(CutMethod,'angle')
            CutCoordList=[CutCoordList;AngleCutPoints{CutList(j,1),CutList(j,2)}];
        else
            error('%s: Cut method unspecified, use either ''distance'' or ''angle''\n',mfilename);
        end
    end
end
    %% Perform cutting: use vision toolbox for optimal performance?
    %using vision toolbox:
%     LineInserter=video.ShapeInserter('Shape','Lines','BorderColor','Custom','CustomBorderColor',1,'Antialiasing',false);
%     CutMask=zeros(ImSize);
%     for j=1:size(CutCoordList,1)
%         CurrentLine=CutCoordList(j,:);
%         CutMask(CurrentLine(1),CurrentLine(2))=1;
%         CutMask(CurrentLine(3),CurrentLine(4))=1;
%         CutMask=step(LineInserter,CutMask,CurrentLine);
%     end

    %using own function
    CutMask=zeros(ImSize);
    for j=1:size(CutCoordList,1)
        CutMask=SimpleLine(CutMask,CutCoordList(j,1:2),CutCoordList(j,3:4),1);
    end
    
    CutMask=CutMask>0;%make sure is logical
    LineStrel=strel('disk',2,0);%maybe use radius 2 for vison toolbox
    CutMask=imdilate(CutMask,getnhood(LineStrel));%thicken cut line to avoid few pixels at end of line still connect objects
end