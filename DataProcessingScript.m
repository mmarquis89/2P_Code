
%%% CREATE ANATOMY STACKS ------------------------------------------------------------------------

fileStr = '*stack_*';

dirPath = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06';
CreateAnatomyStack(dirPath, fileStr, []);

%%% PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05_exp2';
preRegRoutine_MM(parentDir);
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06';
preRegRoutine_MM(parentDir);

%%% REGISTRATION -----------------------------------------------------------------------------------
refVol = 6;
refTrial = 1;

regErrorLog = {};
try
    matlabImReg_2P_session_MM('D:\Dropbox (HMS)\2P Data\Imaging Data\2017_07_28\sid_0',refVol,refTrial);
catch
    regErrorLog{end+1} = '7-28 sid 0';
end

%%% MAKE VIDEOS ------------------------------------------------------------------------------------
frameRate = 25;

vidErrLog = {};
parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_07_28';
try
    makeVids_2P(2, parentDir, frameRate)
catch
    vidErrLog{end+1} = '7-28 sid 2';
end

% Zip raw video frames
zipRawFrames(0, 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06')
zipRawFrames(1, 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06')

% Concatenate
concatErrLog = {};
parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_07_30\_Movies';
try
    concatVids_2P(0, parentDir, frameRate)
catch
    concatErrLog{end+1} = '7-30 sid 0';
end



% COUNT VIDEO FRAMES
% parentDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2017_09_06\_Movies';
% sid = 0;
% vidFiles = dir(fullfile(parentDir, ['sid_', num2str(sid), '*tid*.mp4']));
% vidNames = sort({vidFiles.name});
% nTrials = length(vidNames);
% frameCounts = []; 
% 
% for iTrial = 1:nTrials
%     
%     disp(num2str(iTrial))
%     frameCounts(iTrial) = countFrames(fullfile(parentDir, vidNames{iTrial}));
%     
% end
% 
% % Save frame count log
% save(fullfile(parentDir, ['sid_', num2str(sid) '_frameCountLog.mat']), 'frameCounts')


