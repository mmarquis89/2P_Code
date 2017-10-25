
expDate = '2017_10_22';

%%
%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
fileStr = '*stack2*';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
tic
create_anatomy_stack(parentDir, fileStr, 'Stack2');
disp(['Creating the anatomy stacks took ', round(num2str(toc)) ' sec']);

%%% Archive raw anatomy stacks
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
archiveName = 'AnatomyStacks';
filterString = '*stack*';
system7zip(parentDir, archiveName, '7z', filterString, 1)

% %%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% tic
% preReg_routine_MM(parentDir, 1, expDate);
% disp(['Pre-registration processing took ', round(num2str(toc)) ' sec']);
% 
%%% Archive raw imaging data files
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% archiveName = 'TrialData_sid_1';
% filterString = '*sid_1_t*';
% system7zip(parentDir, archiveName, '7z', filterString, 1)

% 
% %%% REGISTRATION -----------------------------------------------------------------------------------
% refVol = 6;
% refTrial = 10;
% 
% tic
% fileName  = 'sid_0_Chan_1_sessionFile';
% matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
% disp(['Registration took ', round(num2str(toc)) ' sec']);
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
% FRAME_RATE = 25;
% sid = 1;
% 
% %%% Make vids
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
% tic
% make_vids_2P(sid, parentDir, FRAME_RATE)
% disp(['Creating behavior vids took ', round(num2str(toc)) ' sec']);
% %%%
% 
% %%% Archive raw video frames
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
% archiveName = ['sid_', num2str(sid), '_RawFrames'];
% system7zip(parentDir, archiveName, '7z', ['*sid_', num2str(sid), '1_t*'], 1);
% %%%
% 
% %%% Concatenate vids
% tic
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
% concat_vids_2P(sid, parentDir, FRAME_RATE)
% disp(['Concatenating behavior vids took ', round(num2str(toc)) ' sec']);
% %%%


%% =============================================================================================================================


expDate = '2017_10_23';

%%% MAKE VIDEOS ------------------------------------------------------------------------------------
FRAME_RATE = 25;
sid = 0;

%%% Make vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
tic
make_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Creating behavior vids took ', round(num2str(toc)) ' sec']);
%%%
% 
% %%% Archive raw video frames
% parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
% archiveName = ['sid_', num2str(sid), '_RawFrames'];
% system7zip(parentDir, archiveName, '7z', ['*sid_', num2str(sid), '_t*'], 1);
% %%%

%%% Concatenate vids
tic
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
concat_vids_2P(sid, parentDir, FRAME_RATE)
disp(['Concatenating behavior vids took ', round(num2str(toc)) ' sec']);
%%%

%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% sid = 0;
% 
% tic
% preReg_routine_MM(parentDir, sid, expDate);
% disp(['Pre-registration processing took ', round(num2str(toc)) ' sec']);

%%% Archive raw imaging data files
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
archiveName = ['TrialData_sid_', num2str(sid)];
filterString = ['*sid_', num2str(sid), '_t*'];
system7zip(parentDir, archiveName, '7z', filterString, 1)


%%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 6;
refTrial = 10;

tic
fileName  = 'sid_0_Chan_1_sessionFile';
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);

tic
fileName  = 'sid_0_Chan_2_sessionFile';
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);

tic
fileName  = 'sid_0_ChanRatio_sessionFile';
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],fileName,refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);



%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
FRAME_RATE = 25;
trialDuration = 16;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'];
sid = 0;
annotationFileName = [expDate, '_sid_0_Annotation.txt'];

tic
test = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', round(num2str(toc)) ' sec']);