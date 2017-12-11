function normcorre_registration(parentDir, fileName)
%========================================================================================================================================
%
% Uses the NoRMCorre algorithm to apply a rigid motion correction to pre-processed imaging data. It must be saved as a .mat file 
% containing the variables 'trialType', 'origFileNames', 'expDate', and either 'wholeSession' or 'regProduct' (the latter in case I want
% to register files that have already been processed by my old registration algorithm). The imaging data array (which must have 
% dimensions [y, x, plane, volume, trial])is the only variable that is actually used, the rest are just re-saved as a new .mat file in 
% the same directory as the source file.
% 
%  INPUT:
%       parentDir   = the full path to the directory containing the data to be registered
%       fileName    = the name of the input .mat file (this will also affect the output file name)
%
%=========================================================================================================================================

load(fullfile(parentDir, fileName));
if exist('regProduct', 'var')
    wholeSession = single(regProduct);
end
reshapeSize = [size(squeeze(wholeSession(:,:,:,1,1))), size(wholeSession, 4) * size(wholeSession, 5)];
concatSession = reshape(wholeSession, reshapeSize);

% Set parameters
options_rigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
                    'max_shift', [25, 25, 2], ...
                    'init_batch', 100, ...
                    'us_fac', 50 ...
                    ); 
% options_nonRigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
%                     'max_shift', [20, 20, 2], ...
%                     'init_batch', 100, ...
%                     'grid_size', [100, 100] ...
%                     );
                
% Rigid registration
tic; [M, ~, ~, ~] = normcorre(concatSession, options_rigid); toc
regProduct = reshape(M, size(wholeSession));
savefast(fullfile(parentDir, ['rigid_', fileName]), 'regProduct', 'trialType', 'origFileNames', 'expDate', 'scanImageInfo');

% Create and save reference images from registered data
refImages = [];
for iPlane = 1:size(regProduct, 3)
    refImages{iPlane} = squeeze(mean(mean( regProduct(:,:,iPlane,:,:), 4), 5)); % --> [y, x]
end
save(fullfile(sessionDir, ['refImages_', fileName]), 'refImages');


end