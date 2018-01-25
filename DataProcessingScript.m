%==================================================================================================
%%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
%==================================================================================================

expDate = '2018_01_17_exp_1'

fileStr = '*Stack*';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
tic
create_anatomy_stack(parentDir, fileStr, 'Stack');
disp(['Creating the anatomy stacks took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s anatomy stack created in %s min', expDate, num2str(round(toc/60, 1))));
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
%==================================================================================================
% % <3 you Michael!

% -------------------------------------------------------------------------------------------------

expDate = '2018_01_17_exp_2'

parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
tic
sids = 1;

preReg_routine_MM(parentDir, sids, expDate);
disp(['Pre-registration processing took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s pre-processing completed in %s min', expDate, num2str(round(toc/60, 1))));

% -------------------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% NoRMCorre REGISTRATION --------------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDates =  {'2018_01_17_exp_2', ...
            };       
sids = [1 ...
        ];

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp};
    sid = sids(iExp);
    
    parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
    
    fileName = ['sid_', num2str(sid), '_sessionFile.mat'];
    tic; normcorre_registration(parentDir, fileName);
    writeToLog(sprintf('%s_%s NoRMCorre registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
    
    
%     fileName = ['sid_', num2str(sid), '_Chan_1_sessionFile.mat'];
%     tic; normcorre_registration(parentDir, fileName);
%     writeToLog(sprintf('%s_%s NoRMCorre registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
    
    
%     fileName = ['sid_', num2str(sid), '_Chan_2_sessionFile.mat'];
%     tic; normcorre_registration(parentDir, fileName);
%     writeToLog(sprintf('%s_%s NoRMCorre registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
    
end

% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% CREATE BEHAVIOR VIDEOS -------------------------------------------------------------------------
%==================================================================================================

expDate = '2018_01_17_exp_1'

FRAME_RATE = 25;
sid = 1;

%%% Make vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
tic
make_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Creating behavior vids took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s behavior vids created in %s min', expDate, num2str(round(toc/60, 1))));

%%% Concatenate vids
tic
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
concat_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Concatenating behavior vids took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s behavior vids concatenated in %s min', expDate, num2str(round(toc/60, 1))));

%%% Define ROIs for optic flow combined vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
select_video_ROIs(parentDir, sid);
writeToLog(sprintf('%s optic flow ROIs defined', expDate));

%%% Make optic flow combined vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
tic
make_optic_flow_vid(sid, parentDir, FRAME_RATE)
disp(['Creating combined vids took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s combined vids created in %s min', expDate, num2str(round(toc/60, 1))));

%% -------------------------------------------------------------------------------------------------

expDate = '2018_01_17_exp_2'

FRAME_RATE = 25;
sid = 1;

%%% Make vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
tic
make_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Creating behavior vids took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s behavior vids created in %s min', expDate, num2str(round(toc/60, 1))));

%%% Concatenate vids
tic
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
concat_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Concatenating behavior vids took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s behavior vids concatenated in %s min', expDate, num2str(round(toc/60, 1))));

%% % Make optic flow combined vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
tic
make_optic_flow_vid(sid, parentDir, FRAME_RATE)
disp(['Creating combined vids took ', num2str(round(toc/60, 1)) ' min']);
writeToLog(sprintf('%s combined vids created in %s min', expDate, num2str(round(toc/60, 1))));

% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% ARCHIVE FILES-----------------------------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDate = '2017_11_27_exp_1';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 1;

% -------------------------------------------------------------------------------------------------

expDate = '2017_11_27_exp_2';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 2;

%%% Archive raw anatomy stacks
archiveName = 'AnatomyStacks';
filterString = '*Stack*';
system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% Archive raw imaging data files
archiveName = ['TrialData_sid_', num2str(sid)];
filterString = ['*sid_', num2str(sid), '_t*'];
system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% Archive raw video frames
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
archiveName = ['sid_', num2str(sid), '_RawFrames'];
system7zip(parentDir, archiveName, '7z', ['*sid_', num2str(sid), '_t*'], 1);

% -------------------------------------------------------------------------------------------------
%==================================================================================================
%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
%==================================================================================================

FRAME_RATE = 25;
trialDuration = 20;
expDate = '2017_11_29_exp_2';
sid = 0;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
annotationFileName = [expDate, '_Annotation.txt'];

tic
test = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', num2str(round(toc/60, 1)) ' min']);

% -------------------------------------------------------------------------------------------------

%% ==================================================================================================
% REGISTRATION -----------------------------------------------------------------------------------
% preview_trial_movie(wholeSession, 5, 5, [], [], []);
%==================================================================================================
% 
% % -------------------------------------------------------------------------------------------------
% 
% expDate = '2017_11_27_exp_1'
% refVol = 20;
% refTrial = 5;
% sid = 1;
% 
% fileName = ['sid_', num2str(sid), '_Chan_2_sessionFile'];
% % fileName  = ['sid_', num2str(sid), '_sessionFile'];
% tic; matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)],fileName,refVol,refTrial);
% disp(['Registration took ', num2str(round(toc/60, 1)) ' min']);
% writeToLog(sprintf('%s_%s Mattes mutual info registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
% 
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], fileName, [fileName, '_registration_transforms']); 
% newfileName = ['sid_', num2str(sid), '_Chan_1_sessionFile'];
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 
% newfileName = ['sid_', num2str(sid), '_ChanRatio_sessionFile'];
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 

% -------------------------------------------------------------------------------------------------
% 
% expDate = '2017_11_27_exp_2'
% refVol = 19;
% refTrial = 5;
% sid = 2;
% 
% fileName = ['sid_', num2str(sid), '_Chan_2_sessionFile'];
% % fileName  = ['sid_', num2str(sid), '_sessionFile'];
% tic; matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)],fileName,refVol,refTrial);
% disp(['Registration took ', num2str(round(toc/60, 1)) ' min']);
% writeToLog(sprintf('%s_%s Mattes mutual info registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
% 
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], fileName, [fileName, '_registration_transforms']); 
% newfileName = ['sid_', num2str(sid), '_Chan_1_sessionFile'];
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 
% newfileName = ['sid_', num2str(sid), '_ChanRatio_sessionFile'];
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 
% 
% % -------------------------------------------------------------------------------------------------
% 
% expDate = '2017_11_27_exp_3'
% refVol = 0;
% refTrial = 11;
% sid = 1;
% 
% fileName = ['sid_', num2str(sid), '_Chan_2_sessionFile'];
% % fileName  = ['sid_', num2str(sid), '_sessionFile'];
% tic; matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)],fileName,refVol,refTrial);
% disp(['Registration took ', num2str(round(toc/60, 1)) ' min']);
% writeToLog(sprintf('%s_%s Mattes mutual info registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
% 
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], fileName, [fileName, '_registration_transforms']); 
% newfileName = ['sid_', num2str(sid), '_Chan_1_sessionFile'];
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 
% newfileName = ['sid_', num2str(sid), '_ChanRatio_sessionFile'];
% registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 



% %% %%%%%% ATTEMPT TO REGISTER ANATOMY STACKS %%%%%%%%%
% clear all
% dirPath = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_11_27_exp_2';
% fileStr = '*stack_0*.tif';
% outputFilePrefix = [];
% 
% writeToLog(sprintf('Attempting to register 2017_11_27_exp_2 anatomy stacks'));
% tic
% 
% stacks = dir(fullfile(dirPath, fileStr));
% nStacks = 30%length(stacks);
% 
% firstStack = uint16(read_tif(fullfile(dirPath, stacks(1).name)));           % --> [y, x, plane, stackNum, channel]
% currStack = uint16(zeros(size(firstStack)));
% fs = firstStack(:,:,:,1,1);
% allStacks = uint16(zeros([size(fs), nStacks]));                             % --> [y, x, plane, stackNum]
% for iStack = 1:nStacks
%     disp(['Reading ' stacks(iStack).name, '...']);
%     currStack = uint16(read_tif(fullfile(dirPath, stacks(iStack).name)));   % --> [y, x, plane, stackNum, channel]   
%     allStacks(:,:,:, iStack) = currStack(:,:,:,:, 2);                       % --> [y, x, plane, stackNum]
% end
% 
% writeToLog(sprintf('Anatomy stacks loaded in %s min', num2str(round(toc/60, 1))));
% 
% % Set parameters
% options_rigid = NoRMCorreSetParms('d1', size(allStacks, 1), 'd2', size(allStacks, 2), 'd3', size(allStacks, 3), ...
%                     'max_shift', [25, 25, 25], ...
%                     'init_batch', 3 ...
%                     ); 
% % options_nonRigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
% %                     'max_shift', [25, 25, 25], ...
% %                     'init_batch', 3, ...
% %                     'grid_size', [100, 100] ...
% %                     );
% 
% % Rigid
% tic; [M, ~, ~, ~] = normcorre(allStacks, options_rigid); toc
% regProduct = mean(M, 4); % --> [y, x, plane]
% savefast(fullfile(dirPath, ['MeanStack_Rigid']), 'regProduct');
% writeToLog(sprintf('2017_11_27_exp_2 anatomy stack registration complete in %s min', num2str(round(toc/60, 1))));
% 
% % 
% % % Non-rigid
% % tic; [M, ~, ~, ~] = normcorre(concatSession, options_nonRigid); toc
% % regProduct = mean(M, 4); % --> [y, x, plane]
% % savefast(fullfile(dirPath, ['MeanStack_NonRigid']), 'regProduct');   
% %%
% clear all
% dirPath = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_11_27_exp_3';
% fileStr = '*stack_0*.tif';
% outputFilePrefix = [];
% 
% writeToLog(sprintf('Attempting to register 2017_11_27_exp_3 anatomy stacks'));
% 
% stacks = dir(fullfile(dirPath, fileStr));
% nStacks = 10%length(stacks);
% 
% firstStack = uint16(read_tif(fullfile(dirPath, stacks(1).name)));           % --> [y, x, plane, stackNum, channel]
% currStack = uint16(zeros(size(firstStack)));
% fs = firstStack(:,:,:,1,1);
% allStacks = uint16(zeros([size(fs), nStacks]));                             % --> [y, x, plane, stackNum]
% for iStack = 1:nStacks
%     disp(['Reading ' stacks(iStack).name, '...']);
%     currStack = uint16(read_tif(fullfile(dirPath, stacks(iStack).name)));   % --> [y, x, plane, stackNum, channel]   
%     allStacks(:,:,:, iStack) = currStack(:,:,:,:, 2);                       % --> [y, x, plane, stackNum]
% end
% writeToLog(sprintf('Anatomy stacks loaded in %s min', num2str(round(toc/60, 1))));
% 
% % Set parameters
% options_rigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
%                     'max_shift', [25, 25, 25], ...
%                     'init_batch', 2 ...
%                     ); 
% options_nonRigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
%                     'max_shift', [25, 25, 25], ...
%                     'init_batch', 3, ...
%                     'grid_size', [100, 100] ...
%                     );
% 
% % Rigid
% tic; [M, ~, ~, ~] = normcorre(concatSession, options_rigid); toc
% regProduct = mean(M, 4); % --> [y, x, plane]
% savefast(fullfile(dirPath, ['MeanStack_Rigid']), 'regProduct');
% writeToLog(sprintf('2017_11_27_exp_2 anatomy stack registration complete in %s min', num2str(round(toc/60, 1))));
% 
% % % Non-rigid
% % tic; [M, ~, ~, ~] = normcorre(concatSession, options_nonRigid); toc
% % regProduct = mean(M, 4); % --> [y, x, plane]
% % savefast(fullfile(dirPath, ['MeanStack_NonRigid']), 'regProduct');   







