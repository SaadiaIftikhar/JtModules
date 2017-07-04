function lastdir=getlastdir(strRootPath)
% getlastdir retreives the last directory or file name from a path string

    if nargin==0
        strRootPath = {'Z:\Data\Users\YF_DG\20080309165648_M1_080309_YF_DG_batch1_CP004-1de\BATCH\Measurements_Image_FileNames.mat'};
    end

    if ischar(strRootPath)
        lastdir = doit(strRootPath);
    elseif iscell(strRootPath)
        lastdir = cellfun(@doit,strRootPath,'UniformOutput',0);
    else
        error('unknown input type for getlastdir')
    end

end

function strOutput = doit(strInput)

    if strncmp(strInput,'\\',2)
        strInput = strInput(3:end);
    end

    strInput = strrep(strInput,strcat(filesep,filesep),filesep);
    strOutput = strInput;
    
    matFilesepIndices = strfind(strInput, filesep);
    
    if isempty(matFilesepIndices) && strcmp(filesep,'\')
        matFilesepIndices = strfind(strInput, '/');        
    elseif isempty(matFilesepIndices) && strcmp(filesep,'/')
        matFilesepIndices = strfind(strInput, '\');
    end
        
    intPathLength = length(strInput);
    if ~isempty(matFilesepIndices)    
        if matFilesepIndices(end) == intPathLength
            strOutput = strInput(matFilesepIndices(end-1)+1:end-1);
        else
            strOutput = strInput(matFilesepIndices(end)+1:end);  
        end
    end
end