
%%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
fileStr = '*stack_*';
dirPath = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06';
CreateAnatomyStack(dirPath, fileStr, []);

%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06';
preRegRoutine_MM(parentDir);

%%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 6;
refTrial = 1;

matlabImReg_2P_session_MM('D:\Dropbox (HMS)\2P Data\Imaging Data\2017_07_28\sid_0',refVol,refTrial);

%%% MAKE VIDEOS ------------------------------------------------------------------------------------
frameRate = 25;

parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06';
makeVids_2P(0, parentDir, frameRate)

%%% Zip raw video frames
zipRawFrames(0, 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_07_30')

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






