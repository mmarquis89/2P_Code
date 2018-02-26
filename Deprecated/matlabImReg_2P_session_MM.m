function tForms = matlabImReg_2P_session_MM(path, fileName, refVol, refTrial)
%===================================================================================================
% Align all data from imaging session to a single volume (rigid - translation + rotation)
%
% INPUTS:
%       path        = session folder containing unregistered session data file.
%       fileName    = name (no file extension) of the .mat file containing the imaging data/metadata
%       refVol      = the volume number to use as a reference for registration.
%       refTrial    = trial # to pull the reference volume from.
%
% OUTPUTS:
%       tForms      = nTrials x nVols cell array containing the 3dAffine objects specifying the 
%                     registration transformation for each volume.
% 
% Code from AKM, modified by MM
% Last modified: 31-Oct-2017
%===================================================================================================

% Load the unregistered session file
disp('Loading file...')
in = load(fullfile(path, [fileName, '.mat']));
disp('Registering...')

movingFile = in.wholeSession; % [x y plane vol trial]

% Pull out the reference volume
refVolData = squeeze(in.wholeSession(:,:,:,refVol,refTrial)); % --> [y, x, plane]

% Set up registration parameters
[optimizer, ~] = imregconfig('monomodal');
metric = registration.metric.MattesMutualInformation;
metric.UseAllPixels = 1;
% metric.NumberOfSpatialSamples = 50000;
optimizer.GradientMagnitudeTolerance = 1e-4; %1e-3
optimizer.MinimumStepLength = 1e-5;
optimizer.MaximumStepLength = 5e-2;% 0.0625
optimizer.MaximumIterations = 100;
optimizer.RelaxationFactor = 0.5;

% Run registration
tic
tForms = [];
for iTrial=1:size(movingFile,5)
    for iVol=1:size(movingFile,4)
        
        % Isolate current volume and register to reference volume
        currVol = squeeze(movingFile(:,:,:,iVol,iTrial));    % --> [y, x, plane, volume, trial]
        disp(['Volume #', num2str(iVol), ' of ', num2str(size(movingFile,4)), ', Trial #', num2str(iTrial), ' of ' num2str(size(movingFile,5))])
        tForms{iTrial, iVol} = imregtform(currVol, refVolData, 'translation', optimizer, metric, 'PyramidLevels', 2);
        
    end
end
tE_sec=toc;

% Save output
save(fullfile(path, [fileName, '_registration_transforms']),'tForms', 'tE_sec');

end








