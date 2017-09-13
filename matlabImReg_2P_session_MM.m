function matlabImReg_2P_session_MM(path,refVol,refTrial)
%===================================================================================================
% Align all data from imaging session to a single volume (translation)

% Input:
%   path = session folder containeing unregistered cdata files 
%   refVol = the volume number to use as a reference for registration
%   refTrial = trial # to pull the reference volume from
% Output:
%   filename_reg1 (regProduct')
% 
% Code from AKM, modified by MM
% Last modified: 14-Aug-2017
%===================================================================================================
cd(path);

% Load the unregistered session file
files=dir('./*sessionFile.mat');
in=load(files(1).name);

movingFile=in.wholeSession; % [x y plane vol trial]
trialType=in.trialType;  
origFileNames=in.origFileNames;

% refVol is a volume from the baseline period 
% of the first odor trial.
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

% Save output
save('sessionOutfile_Reg1','regProduct','trialType','origFileNames','tE_sec');
    
end