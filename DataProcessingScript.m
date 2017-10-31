
expDate = '2017_10_26_TH_Gal4_exp_2';

%%
% %% CREATE ANATOMY STACKS --------------------------------------------------------------------------
% fileStr = '*Stack*';
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% tic
% create_anatomy_stack(parentDir, fileStr, 'Stack');
% disp(['Creating the anatomy stacks took ', round(num2str(toc)) ' sec']);

% %%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% tic
% sids = 1;
% preReg_routine_MM(parentDir, sids, expDate);
% disp(['Pre-registration processing took ', round(num2str(toc)) ' sec']);

% %%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 8;
refTrial = 10;
sid = 1;
% 
tic
fileName  = ['sid_', num2str(sid), '_sessionFile'];
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)],fileName,refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);
% 
% tic
% fileName  = 'sid_0_Chan_2_sessionFile';
% matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
% disp(['Registration took ', round(num2str(toc)) ' sec']);
% 
% tic
% fileName  = 'sid_0_ChanRatio_sessionFile';
% matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
% disp(['Registration took ', round(num2str(toc)) ' sec']);


%%% MAKE VIDEOS ------------------------------------------------------------------------------------
FRAME_RATE = 25;
sid = 0;

%%% Make vids
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
% tic
% make_vids_2P(sid, parentDir, FRAME_RATE)
% disp(['Creating behavior vids took ', round(num2str(toc)) ' sec']);
%%%

%%% Concatenate vids
% tic
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
% concat_vids_2P(sid, parentDir, FRAME_RATE)
% disp(['Concatenating behavior vids took ', round(num2str(toc)) ' sec']);
% %%%


%%% ARCHIVE FILES ----------------------------------------------------------------------------------
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
sid = 1;

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
%%%




%% =============================================================================================================================




expDate = '2017_10_26_TH_Gal4_exp_2';

%%%
% %% CREATE ANATOMY STACKS --------------------------------------------------------------------------
% fileStr = '*Stack*';
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% tic
% create_anatomy_stack(parentDir, fileStr, 'Stack');
% disp(['Creating the anatomy stacks took ', round(num2str(toc)) ' sec']);

% %%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% tic
% sids = 1;
% preReg_routine_MM(parentDir, sids, expDate);
% disp(['Pre-registration processing took ', round(num2str(toc)) ' sec']);

% %%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 8;
refTrial = 10;
sid = 1;
% 
tic
fileName  = ['sid_', num2str(sid), '_sessionFile'];
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)],fileName,refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);
% 
% tic
% fileName  = 'sid_0_Chan_2_sessionFile';
% matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
% disp(['Registration took ', round(num2str(toc)) ' sec']);
% 
% tic
% fileName  = 'sid_0_ChanRatio_sessionFile';
% matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
% disp(['Registration took ', round(num2str(toc)) ' sec']);


%%% MAKE VIDEOS ------------------------------------------------------------------------------------
FRAME_RATE = 25;
sid = 0;

%%% Make vids
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
% tic
% make_vids_2P(sid, parentDir, FRAME_RATE)
% disp(['Creating behavior vids took ', round(num2str(toc)) ' sec']);
%%%

%%% Concatenate vids
% tic
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
% concat_vids_2P(sid, parentDir, FRAME_RATE)
% disp(['Concatenating behavior vids took ', round(num2str(toc)) ' sec']);
% %%%


%%% ARCHIVE FILES ----------------------------------------------------------------------------------
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
%%%





%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
FRAME_RATE = 25;
trialDuration = 20;
expDate = '2017_10_23';

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'];
sid = 0;
annotationFileName = [expDate, '_sid_0_Annotation.txt'];

tic
test = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', round(num2str(toc)) ' sec']);