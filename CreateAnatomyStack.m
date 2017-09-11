function CreateAnatomyStack(dirPath, fileStr, zipFileName)
%===================================================================================================
% MAKE AVERAGED ANATOMY STACK FROM INDIVIDUAL 2P VOLUMES
% Will average together all stacks in the directory that meet the requirement specified by 'fileStr'
% and save as .tif files the mean stack as well as a max Z-projection. Then, moves all the original 
% stack files to a zipped folder.
%    dirPath = folder with stack files, e.g. 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05'
%    fileStr = filter string for dir() selection to get anatomy stack files, e.g. '*stack_*'. Be 
%              careful to ensure that only the desired files will meet this specification.
%    zipFileName = name to append to 'ZippedAnatomyStacks' for zipfile. Can be [].
%===================================================================================================

stacks = dir(fullfile(dirPath, fileStr));
nStacks = length(stacks);

% Throw error if the files about to be created already exist
assert(exist(fullfile(dirPath, 'MeanMaxZ.tif'), 'file')==0 ...
    && exist(fullfile(dirPath, 'MeanStack.tif'), 'file')==0 ...
    , 'Error: output .tif files already exist in this directory');

% Get sum of all stacks
summedStacks = [];
for iVol = 1:nStacks
    disp(['Reading ' stacks(iVol).name, '...']);
    if iVol == 1
        firstStack = uint32(readTif(fullfile(dirPath,stacks(iVol).name)));
        summedStacks = uint32(firstStack);
    else
        summedStacks = summedStacks + uint32(readTif(fullfile(dirPath,stacks(iVol).name)));
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
imwrite(maxZScaled, fullfile(dirPath, 'MeanMaxZ.tif'));
imwrite(avgStackScaled(:,:,1), fullfile(dirPath, 'MeanStack.tif'));
for iPlane = 2:size(avgStackScaled, 3)
    imwrite(avgStackScaled(:,:,iPlane), fullfile(dirPath, 'MeanStack.tif'), 'writemode','append');
end

% Zip raw anatomy stack images
tic
disp('Mean stack creation complete')
disp('Zipping individual anatomy stacks...')
zipPaths = strcat([dirPath, '\'], {stacks.name});
if isempty(zipFileName)
    zip(fullfile(dirPath,'ZippedAnatomyStacks'), zipPaths);
else
    zip(fullfile(dirPath,['ZippedAnatomyStacks_', zipFileName]), zipPaths);
end
disp(['Zipping complete. Duration = ', num2str(toc), ' sec'])

% Delete raw anatomy stack images if zipping was successful
disp('Deleting raw anatomy stacks...')
delFiles = stacks;
zipDir = dir(fullfile(dirPath, ['ZippedAnatomyStacks_', fileStr(1:end-1), '.zip']));
assert(isempty(zipDir), 'Error: no zipped folder was found for these files');
for iFile = 1:length(delFiles)
    delete(fullfile(dirPath, delFiles(iFile).name));
end
disp('Raw stacks deleted')

end






