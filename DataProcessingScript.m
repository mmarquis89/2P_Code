
expDate = '2017_09_24';

%%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
% fileStr = '*stack1_*';
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% CreateAnatomyStack(parentDir, fileStr, 'Stack');

% fileStr = '*stack2_*';
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% CreateAnatomyStack(parentDir, fileStr, 'Stack2');

%%% Archive raw anatomy stacks
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% archiveName = 'AnatomyStacks';
% filterString = '*stack_*';
% system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
preRegRoutine_MM(parentDir, []);

%%% Archive raw imaging data files
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
archiveName = 'TrialData_sid_0';
filterString = '*sid_0_t*';
system7zip(parentDir, archiveName, '7z', filterString, 1)
% 
% parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
% archiveName = 'TrialData_sid_1';
% filterString = '*sid_1_t*';
% system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 6;
refTrial = 1;
tic
disp('Registering...')
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],refVol,refTrial);
disp(num2str(toc))

%%% MAKE VIDEOS ------------------------------------------------------------------------------------
frameRate = 25;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
makeVids_2P(0, parentDir, frameRate)

%%% Archive raw video frames
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
archiveName = 'sid_0_RawFrames';
system7zip(parentDir, archiveName, '7z', '*sid_0_t*', 1);

%%% Concatenate vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
concatVids_2P(0, parentDir, frameRate)

%%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
frameRate = 25;
trialDuration = 15;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'];
sid = 0;
annotationFileName = [expDate, '_sid_0_Annotation.txt'];

test = processAnvilAnnotations(sid, parentDir, saveDir, annotationFileName, frameRate, trialDuration);






