
%%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
fileStr = '*stack_*';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06';
CreateAnatomyStack(parentDir, fileStr, []);

%%% Archive raw anatomy stacks
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_07_30\ZippedAnatomyStacks_stack';
archiveName = 'AnatomyStacks';
filterString = '*stack*';
system7zip(parentDir, archiveName, '7z', filterString)

%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06';
preRegRoutine_MM(parentDir);

%%% Archive raw imaging data files
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_07_30\TrialData_sid_0+5_OldNames';
archiveName = 'testArchive';
filterString = '*sid_0_t*';
system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 6;
refTrial = 1;

matlabImReg_2P_session_MM('D:\Dropbox (HMS)\2P Data\Imaging Data\2017_07_28\sid_0',refVol,refTrial);

%%% MAKE VIDEOS ------------------------------------------------------------------------------------
frameRate = 25;

parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06';
makeVids_2P(0, parentDir, frameRate)

%%% Archive raw video frames
parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06';
archiveName = 'sid_0_RawFrames';
system7zip(parentDir, archiveName, '7z', '*sid_0_t*', 1);

%%% Concatenate vids
parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06\_Movies';
concatVids_2P(0, parentDir, frameRate)

%%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
frameRate = 25;
trialDuration = 10;

parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06\_Movies';
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06\sid_0';
sid = 0;
annotationFileName = '2017_09_06_Annotation_sid_0.txt';

test = processAnvilAnnotations(sid, parentDir, saveDir, annotationFileName, frameRate, trialDuration);






