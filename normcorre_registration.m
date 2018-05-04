function normcorre_registration(parentDir, fileName, varargin)
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
% OPTIONAL NAME-VALUE PAIR INPUTS:
%       
%       'OutputDir' = (default: parentDir) the directory to save the processed data in
%=========================================================================================================================================

addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutputDir', parentDir);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;

wholeSession = load(fullfile(parentDir, fileName));
reshapeSize = [size(squeeze(wholeSession(:,:,:,1,1))), size(wholeSession, 4) * size(wholeSession, 5)];
concatSession = reshape(wholeSession, reshapeSize);

% Make session folder for new files if necessary
if ~isdir(outputDir)
    mkdir(outputDir)
end

% Set parameters
options_rigid = NoRMCorreSetParms('d1', sz(1), 'd2', sz(2), 'd3', sz(3), ...
                    'max_shift', [25, 25, 2], ...
                    'init_batch', 100, ...
                    'us_fac', 50, ...
                    ); 
% options_nonRigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
%                     'max_shift', [20, 20, 2], ...
%                     'init_batch', 100, ...
%                     'grid_size', [100, 100] ...
%                     );
                
% Rigid registration
tic; [M, ~, ~, ~] = normcorre_batch(m, options_rigid); toc
wholeSession = reshape(M, size(wholeSession));

% Save registered data
save(fullfile(parentDir, ['rigid_', fileName]), 'wholeSession')

% Create and save reference images from registered data
refImages = [];
for iPlane = 1:size(regProduct, 3)
    refImages{iPlane} = squeeze(mean(mean( regProduct(:,:,iPlane,:,:), 4), 5)); % --> [y, x]
end
 savefast(fullfile(parentDir, ['refImages_Reg.mat']), 'refImages');


end