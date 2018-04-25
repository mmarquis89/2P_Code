function registration_transform(parentDir, inputDataFileName, transformDataFileName)
%===================================================================================================
% Register imaging data using previously calculated registration transforms
%
% INPUTS:
%       parentDir   = session folder containing unregistered session data file.
%       fileName    = name (no file extension) of the .mat file containing the imaging data/metadata
%       refVol      = the volume number to use as a reference for registration.
%       refTrial    = trial # to pull the reference volume from.
%
% Code from AKM, modified by MM
% Last modified: 31-Oct-2017
%===================================================================================================

% Load the unregistered session file
disp('Loading files...')
inputData = load(fullfile(parentDir, [inputDataFileName, '.mat']));
sessionData = inputData.wholeSession;
trialType = inputData.trialType;  
origFileNames = inputData.origFileNames;
expDate = inputData.expDate;
scanImageInfo = inputData.scanImageInfo;

% Load the transforms file
transformData = load(fullfile(parentDir, [transformDataFileName, '.mat']));
tForms = transformData.tForms;
tE_sec = transformData.tE_sec;

% Initialize empty array
regProduct=uint16(zeros(size(sessionData)));

% Transform the raw session data
for iTrial = 1:size(transformData.tForms, 1)
    for iVol = 1:size(transformData.tForms, 2)
        currVol = squeeze(sessionData(:,:,:, iVol, iTrial));
        regProduct(:,:,:, iVol, iTrial) = imwarp(currVol, tForms{iTrial, iVol}, 'OutputView', imref3d(size(currVol))); 
    end    
end

% Remove edge pixels to crop out any registration artifacts
yrange = 3:size(regProduct,1)-2;
xrange = 3:size(regProduct,2)-2;
regProduct = regProduct(yrange,xrange,:,:,:); % --> [y, x, plane, volume, trial]

% Re-offset so minimum value = 1
regProduct = (regProduct - min(regProduct(:))) + 1;
disp(['nan count = ', num2str(sum(isnan(regProduct (:))))]);
disp(['inf count = ', num2str(sum(isinf(regProduct (:))))]);

% Save output
outputFileName = [inputDataFileName, '_Reg1'];
disp(['Saving registered output as ', outputFileName, '.mat...'])
save(fullfile(parentDir, outputFileName),'regProduct','trialType','origFileNames', 'scanimageInfo', 'tE_sec','expDate');
disp('Saving complete')

end
