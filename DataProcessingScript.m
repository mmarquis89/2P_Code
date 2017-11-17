%==================================================================================================
%%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
%==================================================================================================

expDate = '2017_11_15_exp_3'

fileStr = '*Stack*';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
tic
create_anatomy_stack(parentDir, fileStr, 'Stack');
disp(['Creating the anatomy stacks took ', round(num2str(toc)) ' sec']);

% -------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% CREATE BEHAVIOR VIDEOS -------------------------------------------------------------------------
%==================================================================================================


expDate = '2017_11_15_exp_2'

FRAME_RATE = 25;
sid = 1;

%%% Make vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
tic
make_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Creating behavior vids took ', round(num2str(toc)) ' sec']);

%%% Concatenate vids
tic
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
concat_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Concatenating behavior vids took ', round(num2str(toc)) ' sec']);

% -------------------------------------------------------------------------------------------------


%==================================================================================================
%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
%==================================================================================================
% % <3 you Michael!

expDate = '2017_11_15_exp_2'

parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
tic
sids = 1;

preReg_routine_MM(parentDir, sids, expDate);
disp(['Pre-registration processing took ', round(num2str(toc)) ' sec']);

% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% REGISTRATION -----------------------------------------------------------------------------------
% preview_trial_movie(myData.wholeSession, 12, 11, [], [], []);
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDate = '2017_11_15_exp_2'
refVol = 90;
refTrial = 11;
sid = 1;

tic
fileName = ['sid_', num2str(sid), '_Chan_2_sessionFile'];
% fileName  = ['sid_', num2str(sid), '_sessionFile'];
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)],fileName,refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);

registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], fileName, [fileName, '_registration_transforms']); 
newfileName = ['sid_', num2str(sid), '_Chan_1_sessionFile'];
registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 
newfileName = ['sid_', num2str(sid), '_ChanRatio_sessionFile'];
registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 

% -------------------------------------------------------------------------------------------------

expDate = '2017_11_15_exp_3'
refVol = 31;
refTrial = 21;
sid = 0;

tic
fileName = ['sid_', num2str(sid), '_Chan_2_sessionFile'];
% fileName  = ['sid_', num2str(sid), '_sessionFile'];
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)],fileName,refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);

registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], fileName, [fileName, '_registration_transforms']); 
newfileName = ['sid_', num2str(sid), '_Chan_1_sessionFile'];
registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 
newfileName = ['sid_', num2str(sid), '_ChanRatio_sessionFile'];
registration_transform(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)], newfileName, [fileName, '_registration_transforms']); 

%==================================================================================================
%% ARCHIVE FILES-----------------------------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDate = '2017_11_14_exp_3';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 0;

%%% Archive raw video frames
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
archiveName = ['sid_', num2str(sid), '_RawFrames'];
system7zip(parentDir, archiveName, '7z', ['*sid_', num2str(sid), '_t*'], 1);

% -------------------------------------------------------------------------------------------------

expDate = '2017_11_15';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 3;

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

expDate = '2017_11_15_exp_2';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 1;

%%% Archive raw imaging data files
archiveName = ['TrialData_sid_', num2str(sid)];
filterString = ['*sid_', num2str(sid), '_t*'];
system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% Archive raw video frames
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
archiveName = ['sid_', num2str(sid), '_RawFrames'];
system7zip(parentDir, archiveName, '7z', ['*sid_', num2str(sid), '_t*'], 1);

% -------------------------------------------------------------------------------------------------

expDate = '2017_11_15_exp_3';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 0;

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
trialDuration = 30;
expDate = '2017_11_15_exp_3';
sid = 0;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
annotationFileName = [expDate, '_sid_', num2str(sid), '_Annotation.txt'];

tic
test = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', round(num2str(toc)) ' sec']);

