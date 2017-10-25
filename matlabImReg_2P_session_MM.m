function matlabImReg_2P_session_MM(path, fileName, refVol, refTrial)
%===================================================================================================
% Align all data from imaging session to a single volume (translation)
%
% INPUTS:
%       path = session folder containing unregistered cdata files.
%       fileName = name (no file extension) of the .mat file containing the imaging data/metadata
%       refVol = the volume number to use as a reference for registration.
%       refTrial = trial # to pull the reference volume from.
%
% OUTPUTS:
%       filename_reg1 (regProduct')
% 
% Code from AKM, modified by MM
% Last modified: 23-Oct-2017
%===================================================================================================
% cd(path);

% Load the unregistered session file
disp('Loading file...')
in = load(fullfile(path, [fileName, '.mat']));
disp('Registering...')

movingFile = in.wholeSession; % [x y plane vol trial]
trialType = in.trialType;  
origFileNames = in.origFileNames;
expDate = in.expDate;


% Pull out the reference volume
refVolData=squeeze(in.wholeSession(:,:,:,refVol,refTrial));

% Initialize empty array
regProduct=uint16(zeros(size(movingFile)));
[optimizer, metric]=imregconfig('monomodal');
metric=registration.metric.MattesMutualInformation; % Optimizer won't converge sometimes w/mean squares metric, for reasons I can't understand...
% optimizer.MaximumIterations=100;
% optimizer.GradientMagnitudeTolerance=1e-3;
% optimizer.MaximumStepLength=5e-2;

tic
for iTrial=1:size(movingFile,5)
    for iVol=1:size(movingFile,4)
        
        mVol=squeeze(movingFile(:,:,:,iVol,iTrial));
        disp(['Volume #', num2str(iVol), ' of ', num2str(size(movingFile,4)), ', Trial #', num2str(iTrial), ' of ' num2str(size(movingFile,5))])
        transProd=imregister(mVol,refVolData,'translation',optimizer,metric,'PyramidLevels',2);
        regProduct(:,:,:,iVol,iTrial)=transProd;
        
    end
end
tE_sec=toc;

% Remove edge pixels to crop out any registration artifacts
yrange = 3:size(regProduct,1)-2;
xrange = 3:size(regProduct,2)-2;
regProduct = regProduct(yrange,xrange,:,:);

% Save output
save(fullfile(path, [fileName, '_Reg1']),'regProduct','trialType','origFileNames','tE_sec','expDate');
    
end