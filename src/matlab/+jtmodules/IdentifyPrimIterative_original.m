classdef IdentifyPrimIterative_original < handle

properties (Constant)

    VERSION = '0.0.3'
end

methods (Static)    

    
function [imFinalObjects, figure] =  main(input_image, CuttingPasses, ...
                        ThresholdCorrection, ThresholdRange, SolidityThres, ...
                        FormFactorThres, LowerSizeThres, UpperSizeThres, ...
                        smoothingDiskSize, WindowSize, PerimSegEqSegment, ...
                        PerimSegEqRadius, LowerSizeCutThres, plot)
 
 
ImageName = 'Image';
ModuleName = 'PrimaryIterative';
ObjectName = 'Object'; 
DebugMode = 'Off'; 
TestMode2 = 'false';
TestMode = 'false'; 
 
 handles.Settings.CurrentModule = 'PrimaryIterative'; 
 handles.Settings.CurrentModuleNum = 1;
 handles.Settings.ModuleName = 'PrimaryIterative'; 
 handles.Settings.ObjectName = 'Object'; 
 handles.Settings.ImageName = 'Image';  
 handles.Settings.ThresholdCorrection = ThresholdCorrection; 
 handles.Settings.ThresholdRange = ThresholdRange; 
 handles.Settings.CuttingPasses = CuttingPasses; 
 handles.Settings.DebugMode = 'Off'; 
 handles.Settings.SolidityThres = SolidityThres; 
 handles.Settings.FormFactorThres = FormFactorThres; 
 handles.Settings.UpperSizeThres = UpperSizeThres; 
 handles.Settings.LowerSizeThres = LowerSizeThres; 
 handles.Settings.LowerSizeCutThres = LowerSizeCutThres; 
 handles.Settings.TestMode2 = 'false'; 
 handles.Settings.WindowSize = WindowSize; 
 handles.Settings.smoothingDiskSize = smoothingDiskSize; 
 handles.Settings.PerimSegEqRadius = PerimSegEqRadius; 
 handles.Settings.PerimSegEqSegment = PerimSegEqSegment; 
 handles.Settings.TestMode = 'false'; 
 handles.Settings.input_image = input_image;


%%%%%%%%%%%%%%%%%%%%
% IMAGE ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%

Threshold = 'Otsu Global';
pObject = '10%';

%%% Chooses the first word of the method name (removing 'Global' or 'Adaptive').
ThresholdMethod = strtok(Threshold);
%%% Checks if a custom entry was selected for Threshold, which means we are using an incoming binary image rather than calculating a threshold.
if isempty(strmatch(ThresholdMethod,{'Otsu','MoG','Background','RobustBackground','RidlerCalvard','All','Set'},'exact'))
    %if ~(strncmp(Threshold,'Otsu',4) || strncmp(Threshold,'MoG',3) || strfind(Threshold,'Background') ||strncmp(Threshold,'RidlerCalvard',13) || strcmp(Threshold,'All') || strcmp(Threshold,'Set interactively'))
    if isnan(str2double(Threshold))
        GetThreshold = 0;
        BinaryInputImage = CPretrieveimage(handles,Threshold,ModuleName,'MustBeGray','CheckScale');
    else
        GetThreshold = 1;
    end
else
    GetThreshold = 1;
end

%%% Checks that the Min and Max threshold bounds have valid values
% index = strfind(ThresholdRange,',');
% if isempty(index)
%     error(['Image processing was canceled in the ', ModuleName, ' module because the Min and Max threshold bounds are invalid.'])
% end
% MinimumThreshold = ThresholdRange(1:index-1);
% MaximumThreshold = ThresholdRange(index+1:end);

MinimumThreshold = min(ThresholdRange);
MaximumThreshold = max(ThresholdRange);

if GetThreshold
    input_image(input_image > quantile(input_image(:), 0.999)) = quantile(input_image(:), 0.999);
    [handles,OrigThreshold] = CPthreshold(handles,Threshold,pObject,MinimumThreshold,MaximumThreshold,ThresholdCorrection,input_image,ImageName,ModuleName,ObjectName);
else
    OrigThreshold = 0;
end

%%% Threshold intensity image
ThreshImage = zeros(size(input_image), 'double');
ThreshImage(input_image > OrigThreshold) = 1;

%%% Fill holes in objects
imInputObjects = imfill(double(ThreshImage),'holes');

if ~isempty(imInputObjects)
    
    %-------------------------------------------
    % Select objects in input image for cutting
    %-------------------------------------------
    
    imObjects = zeros([size(imInputObjects),CuttingPasses]);
    imSelected = zeros([size(imInputObjects),CuttingPasses]);
    imCutMask = zeros([size(imInputObjects),CuttingPasses]);
    imCut = zeros([size(imInputObjects),CuttingPasses]);
    imNotCut = zeros([size(imInputObjects),CuttingPasses]);
    objFormFactor = cell(CuttingPasses,1);
    objSolidity = cell(CuttingPasses,1);
    objArea = cell(CuttingPasses,1);
    cellPerimeterProps = cell(CuttingPasses,1);
    
    for i = 1:CuttingPasses
        
        if i==1
            imObjects(:,:,i) = imInputObjects;
        else
            imObjects(:,:,i) = imCut(:,:,i-1);
        end
        
        % Measure basic area/shape features
        props = regionprops(logical(imObjects(:,:,i)),'Area','Solidity','Perimeter');
        
        % Features used for object selection
        objSolidity{i} = cat(1,props.Solidity);
        objArea{i} = cat(1,props.Area);
        tmp = log((4*pi*cat(1,props.Area)) ./ ((cat(1,props.Perimeter)+1).^2))*(-1);%make values positive for easier interpretation of parameter values
        tmp(tmp<0) = 0;
        objFormFactor{i} = tmp;
        
        % Select objects based on these features (user defined thresholds)
        obj2cut = objSolidity{i} < SolidityThres & objFormFactor{i} > FormFactorThres ...
            & objArea{i} > LowerSizeThres & objArea{i} < UpperSizeThres;
        objNot2cut = ~obj2cut;
        
        objSelected = zeros(size(obj2cut));
        objSelected(obj2cut) = 1;
        objSelected(objNot2cut) = 2;
        imSelected(:,:,i) = rplabel(logical(imObjects(:,:,i)),[],objSelected);
        
        % Create mask image with objects selected for cutting
        imObj2Cut = zeros(size(input_image));
        imObj2Cut(imSelected(:,:,i)==1) = 1;
        
        % Store remaining objects that are omitted from cutting
        tmp = zeros(size(input_image));
        tmp(imSelected(:,:,i)==2) = 1;
        imNotCut(:,:,i) = logical(tmp);
        
        
        %-------------
        % Cut objects
        %-------------
        
        % Smooth image
        SmoothDisk = getnhood(strel('disk',double(smoothingDiskSize),0));%minimum that has to be done to avoid problems with bwtraceboundary
        imObj2Cut = bwlabel(imdilate(imerode(imObj2Cut,SmoothDisk),SmoothDisk));
        
        % In rare cases the above smoothing approach creates new, small
        % objects that cause problems. Let's remove them.
        props = regionprops(logical(imObj2Cut),'Area');
        objArea2 = cat(1,props.Area);
        obj2remove = find(objArea2 < LowerSizeThres);
        for j = 1:length(obj2remove)
            imObj2Cut(imObj2Cut==obj2remove(j)) = 0;
        end
        imObj2Cut = bwlabel(imObj2Cut);
        
        % Separate clumped objects along watershed lines
        
        % Note: PerimeterAnalysis cannot handle holes in objects (we may
        % want to implement this in case of big clumps of many objects).
        % Sliding window size is linked to object size. Small object sizes
        % (e.g. in case of images acquired with low magnification) limits
        % maximal size of the sliding window and thus sensitivity of the
        % perimeter analysis.
        
        % Perform perimeter analysis
        cellPerimeterProps{i} = PerimeterAnalysis(imObj2Cut,WindowSize);
        
        % This parameter limits the number of allowed concave regions.
        % It can serve as a safety measure to prevent runtime problems for
        % very complex objects.
        % This could become an input argument in the future!?
        numRegionTheshold = 30;
        
        % Perform the actual segmentation
        if strcmp(DebugMode, 'On')
            imCutMask(:,:,i) = PerimeterWatershedSegmentation(imObj2Cut,input_image,cellPerimeterProps{i},PerimSegEqRadius,PerimSegEqSegment,LowerSizeCutThres, numRegionTheshold, 'debugON');
        else
            imCutMask(:,:,i) = PerimeterWatershedSegmentation(imObj2Cut,input_image,cellPerimeterProps{i},PerimSegEqRadius,PerimSegEqSegment,LowerSizeCutThres, numRegionTheshold);
        end
        imCut(:,:,i) = bwlabel(imObj2Cut.*~imCutMask(:,:,i));
        
        
        %------------------------------
        % Display intermediate results
        %------------------------------
        
        drawnow
        
        % Create overlay images
        imOutlineShapeSeparatedOverlay = input_image;
        B = bwboundaries(imCut(:,:,i));
        imCutShapeObjectsLabel = label2rgb(bwlabel(imCut(:,:,i)),'jet','k','shuffle');
        
        % GUI
%         if CPisHeadless == false
%             tmpSelected = (imSelected(:,:,i));
%             ThisModuleFigureNumber = handles.Current.(['FigureNumberForModule',CurrentModule]);
%             CPfigure(handles,'Image',ThisModuleFigureNumber);
%             subplot(2,2,2); CPimagesc(logical(tmpSelected==1),handles);
%             title(['Cut lines on selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
%             hold on
%             L = bwboundaries(logical(imCutMask(:,:,i)), 'noholes');
%             for l = 1:length(L)
%                 line = L{l};
%                 plot(line(:,2), line(:,1), 'r', 'LineWidth', 3);
%             end
%             hold off
%             freezeColors
%             subplot(2,2,1); CPimagesc(imSelected(:,:,i),handles); colormap('jet');
%             title(['Selected original objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
%             freezeColors
%             subplot(2,2,3); CPimagesc(imOutlineShapeSeparatedOverlay,handles);
%             title(['Outlines of Separated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
%             hold on
%             for k = 1:length(B)
%                 boundary = B{k};
%                 plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
%             end
%             hold off
%             freezeColors
%             subplot(2,2,4); CPimagesc(imCutShapeObjectsLabel,handles);
%             title(['Separated objects, cycle # ',num2str(handles.Current.SetBeingAnalyzed)]);
%             freezeColors
%         end
    end
    
    %-----------------------------------------------
    % Combine objects from different cutting passes
    %-----------------------------------------------
    
    imCut = logical(imCut(:,:,CuttingPasses));
    
    if ~isempty(imCut)
        imErodeMask = bwmorph(imCut,'shrink',inf);
        imDilatedMask = IdentifySecPropagateSubfunction(double(imErodeMask), double(input_image),imCut,1);
    end
    
    imNotCut = logical(sum(imNotCut,3));% Retrieve objects that were not cut
    imFinalObjects = bwlabel(logical(imDilatedMask + imNotCut));
    
else
    
    cellPerimeterProps = {};
    imFinalObjects = zeros(size(imInputObjects));
    imObjects = zeros([size(imInputObjects),CuttingPasses]);
    imSelected = zeros([size(imInputObjects),CuttingPasses]);
    imCutMask = zeros([size(imInputObjects),CuttingPasses]);
    imCut = zeros([size(imInputObjects),CuttingPasses]);
    imNotCut = zeros([size(imInputObjects),CuttingPasses]);
    objFormFactor = cell(CuttingPasses,1);
    objSolidity = cell(CuttingPasses,1);
    objArea = cell(CuttingPasses,1);
    cellPerimeterProps = cell(CuttingPasses,1);
    
end

imFinalObjects = int32(imFinalObjects);

if plot
        plots = { ...
                jtlib.plotting.create_intensity_image_plot(imFinalObjects, 'ul')};
                figure = jtlib.plotting.create_figure(plots);
else
    figure = '';
end

end

end
end