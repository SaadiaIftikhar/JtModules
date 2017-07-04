function FilteredMediafilenames = getFilteredMediafilenamesCP3D(selSpecifiedPathname,strCondition,selChannel,ZPositions)
% This function prefilters image names by a user input string and the
% channel number. Output will be a cell, where rows are different image
% sets and columns are Z-planes of one image set. 
%
% ZPOSITIONS can be either a 1xX vector or 'all' to select all planes or []
% to ignore Z position information and treat all images equally
%
% Note that this is a modified version of a a simplified version of the initializeCP3DStack.

% select subset of files based upon string. note that this allows to easily
% select wells of positive controls (amongst other appliations)
FileList = CPretrievemediafilenames(selSpecifiedPathname,char(strCondition),'', 'Regular','Image')';

% Filter for channel
[intChannel] = cellfun(@(x) check_image_channel(x), FileList,'UniformOutput',false);
intChannel = cell2mat(intChannel);
f = intChannel == selChannel;
FileList = FileList(f);



% in case that no Z Position has been defined (is empty []), do a trival
% case handling, where all names are considered equal and quit function
if isempty(ZPositions)
    FilteredMediafilenames = FileList;
    return;
else
    % otherwise obtain information for sorting according to Z
    [subRCTS, ~, intZPos] = obtainRCTS(FileList);
    % select all Z-planes, if wished
    if ischar(ZPositions)
        if strcmp(ZPositions,'all') == true
            ZPositions=unique(intZPos); % if 'all' is specified, select all z planes in ascending order.
        end
    end
end

% filter for selected Z planes
f = ismember(intZPos,ZPositions);
subRCTS = subRCTS(f,:);
intZPos = intZPos(f);
FileList = FileList(f);

% Create Output

% Trivial solution if only one z-plane is required
if length(ZPositions) == 1
    FilteredMediafilenames = FileList;
% Otherwise group images of one stack according to defined z-plane    
else
    %%% Obtain IDs of images to group in one stack %%%
    [usubRCTS uPos]= unique(subRCTS);       % obtain unique subRCTS
    [~, sortIX]=sort(uPos);                 % restore order of initial file-list > usubRCTSi
    usubRCTSi = usubRCTS(sortIX);
    
    % Check if correct amount of images
    ImgPerRCTS = histc(subRCTS,usubRCTS); % Note that edge vector must be monotonically non-decreasing (thus usubRCTS rather than usubRCTSi is used.
    if any(ImgPerRCTS ~= length(ZPositions)) == true
        error('At least one image set does not have all specified Z-planes');
    end
    clear usubRCTS;  % clear usubRCTS so that it will not be available later, making sure that only unsorted usubRCTSi will be used in later code.
    
    % Initialize Filenames
    dimImageStacks = size(usubRCTSi,1);
    dimImagesPerStack = size(ZPositions,2);
    
    FilteredMediafilenames = cell(dimImageStacks,dimImagesPerStack);
    for l=1:dimImageStacks
        for r = 1: dimImagesPerStack
            currentusubRCTS = usubRCTSi(l,1);
            f = (subRCTS == currentusubRCTS) & (intZPos==ZPositions(1,r));
            FilteredMediafilenames(l,r)=FileList(f,1);
        end
    end
end
end