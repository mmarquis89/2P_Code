function create_anatomy_stack(dirPath, fileStr, outputFilePrefix)
%===================================================================================================
% MAKE AVERAGED ANATOMY STACK FROM INDIVIDUAL 2P VOLUMES
% Will average together all stacks in the directory that meet the requirement specified by 'fileStr'
% and save as .tif files the mean stack as well as a max Z-projection. 
% 
% Inputs:
%    dirPath = folder with stack files, e.g. 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05'
%    fileStr = filter string for dir() selection to get anatomy stack files, e.g. '*stack_*'. Be 
%              careful to ensure that only the desired files will meet this specification.
%    outputFilePrefix = name to prepend to the names of output files. Can be [].
%===================================================================================================

stacks = dir(fullfile(dirPath, fileStr));
nStacks = length(stacks);

% Add underscore to prefix if one is provided
if ~isempty(outputFilePrefix)
   outputFilePrefix = [outputFilePrefix, '_']; 
end

% Throw error if the files about to be created already exist
assert(exist(fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ.tif']), 'file')==0 ... 
    && exist(fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']), 'file')==0 ...
    , 'Error: output .tif files already exist in this directory');

% Get sum of all stacks
summedStacks = [];
for iVol = 1:nStacks
    disp(['Reading ' stacks(iVol).name, '...']);
    if iVol == 1
        firstStack = uint32(read_tif(fullfile(dirPath,stacks(iVol).name)));
        summedStacks = uint32(firstStack);
    else
        summedStacks = summedStacks + uint32(read_tif(fullfile(dirPath,stacks(iVol).name)));
    end
end
% disp(['Max summed value = ', num2str(max(summedStacks(:)))]);

% Convert to double
summedStacks = double(summedStacks);

% Calculate mean by dividing by total number of stacks and scale from 0-1
avgStack = summedStacks ./ nStacks;
avgStackScaled = avgStack ./ max(avgStack(:));

% Get max Z-projection of averaged stack and scale intensity range to 0-1
maxZ = max(avgStack, [], 3);
maxZScaled = maxZ ./ max(maxZ(:));

% Write averaged stack and Z-projection to .tif files after checking for existing files
imwrite(maxZScaled, fullfile(dirPath, [outputFilePrefix, 'MeanMaxZ.tif']));
imwrite(avgStackScaled(:,:,1), fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']));
for iPlane = 2:size(avgStackScaled, 3)
    imwrite(avgStackScaled(:,:,iPlane), fullfile(dirPath, [outputFilePrefix, 'MeanStack.tif']), 'writemode','append');
end







