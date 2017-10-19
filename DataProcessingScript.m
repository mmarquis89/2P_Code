
expDate = '2017_10_17_25A01_testing';


%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
fileStr = '*stack*';
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
tic
create_anatomy_stack(parentDir, fileStr, 'Stack');
disp(['Creating the anatomy stacks took ', round(num2str(toc)) ' sec']);

%%% Archive raw anatomy stacks
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
archiveName = 'AnatomyStacks';
filterString = '*stack*';
system7zip(parentDir, archiveName, '7z', filterString, 1)

%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
tic
preReg_routine_MM(parentDir, [], expDate);
disp(['Pre-registration processing took ', round(num2str(toc)) ' sec']);

%%% Archive raw imaging data files
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
archiveName = 'TrialData_sid_0';
filterString = '*sid_0_t*';
system7zip(parentDir, archiveName, '7z', filterString, 1)


%%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 6;
refTrial = 10;

disp('Registering...')
tic
matlabImReg_2P_session_MM(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'],refVol,refTrial);
disp(['Registration took ', round(num2str(toc)) ' sec']);

%%% MAKE VIDEOS ------------------------------------------------------------------------------------
frameRate = 25;

%%% Make vids
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
tic
make_vids_2P(0, parentDir, frameRate)
disp(['Creating behavior vids took ', round(num2str(toc)) ' sec']);
%%%

%%% Archive raw video frames
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
archiveName = 'sid_0_RawFrames';
system7zip(parentDir, archiveName, '7z', '*sid_0_t*', 1);
%%%

%%% Concatenate vids
tic
parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
concat_vids_2P(0, parentDir, frameRate)
disp(['Concatenating behavior vids took ', round(num2str(toc)) ' sec']);
%%%

%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
frameRate = 25;
trialDuration = 16;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_0'];
sid = 0;
annotationFileName = [expDate, '_sid_0_Annotation.txt'];

tic
test = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, frameRate, trialDuration);
disp(['Processing Anvil annotations took ', round(num2str(toc)) ' sec']);