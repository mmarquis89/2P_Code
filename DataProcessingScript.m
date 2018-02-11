%==================================================================================================
%%% EXTRACT SAMPLE FRAMES FROM ANATOMY STACKS -----------------------------------------------------
%==================================================================================================

expDates = {'2018_02_07_exp_1' ...
            '2018_02_07_exp_2' ...
            '2018_02_07_exp_3' ...
            }

targetPlanes = [ 340 ...
                 155 ...
                 90 ...
    ];

fileStr = '*Stack_0*.tif';

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    targetPlane = targetPlanes(iExp);
    
    dirPath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    disp(['Extracting sample frames from anatomy stacks...']);
    extract_sample_frames(dirPath, fileStr, targetPlane);
    writeToLog(sprintf('%s sample frames extracted in %s min', expDate, num2str(round(toc/60, 1))));
    disp(['Extracting sample frames took ', num2str(round(toc/60, 1)) ' min']);
end
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% % CREATE AVERAGE FLUORESCENCE VIDEOS ------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDates = {'2018_02_07_exp_1' ...
            '2018_02_07_exp_2' ...
            '2018_02_07_exp_3' ...
            }

sids = [ 0 ...
         0 ...
         0 ...
         0 ...
         0 ...
    ];
for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    sid = sids(iExp);
    
    tic
    parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    disp('Creating average fluorescence videos...');
    create_average_fluorescence_vid(parentDir, sid)
    writeToLog(sprintf('%s average fluorescence vids made in %s min', expDate, num2str(round(toc/60, 1))));
    disp(['Creating average fluorescence videos took ', num2str(round(toc/60, 1)) ' min']);
end
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% % CREATE ANATOMY STACKS --------------------------------------------------------------------------
%==================================================================================================

expDates = {'2018_02_07_exp_1' ...
            '2018_02_07_exp_2' ...
            '2018_02_07_exp_3' ...
            }

fileStr = '*Stack*.tif';

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    
    parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    create_anatomy_stack(parentDir, fileStr, 'Stack');
    disp(['Creating the anatomy stacks took ', num2str(round(toc/60, 1)) ' min']);
    writeToLog(sprintf('%s anatomy stack created in %s min', expDate, num2str(round(toc/60, 1))));
end
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% % PRE-REGISTRATION PROCESSING --------------------------------------------------------------------
%==================================================================================================
% % <3 you Michael!

% -------------------------------------------------------------------------------------------------

expDates = {'2018_02_07_exp_1' ...
            '2018_02_07_exp_2' ...
            '2018_02_07_exp_3' ...
            }

sids = [ 0 ...
         0 ...
         0 ...
         0 ...
         0 ...
    ];

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp};
    sid = sids(iExp);
    
    parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    preReg_routine_MM(parentDir, sid, expDate) 
    writeToLog(sprintf('%s_%s pre-processing completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
end

% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% % NoRMCorre REGISTRATION --------------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDates = {'2018_02_07_exp_1' ...
            '2018_02_07_exp_2' ...
            '2018_02_07_exp_3' ...
            }

sids = [ 0 ...
         0 ...
         0 ...
         0 ...
         0 ...
    ];

% fileNames = repmat(['sid_', num2str(sid), '_sessionFile.mat'], 1, 2);
fileNames = { 'sid_0_sessionFile.mat', ...
              'sid_0_sessionFile.mat', ...
              'sid_0_sessionFile.mat', ...
              'sid_0_sessionFile.mat', ...
              'sid_0_sessionFile.mat', ...
            }

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp};
    sid = sids(iExp);
    fileName = fileNames{iExp};
    
    parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
    
    tic; normcorre_registration(parentDir, fileName);
    writeToLog(sprintf('%s_%s NoRMCorre registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
    
end

% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% CREATE BEHAVIOR VIDEOS -------------------------------------------------------------------------
%==================================================================================================

expDates = {'2018_02_07_exp_1' ...
            '2018_02_07_exp_2' ...
            '2018_02_07_exp_3' ...
            }

sids = [ 0 ...
         0 ...
         0 ...
         0 ...
         0 ...
    ];

FRAME_RATE = 25;

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp};
    sid = sids(iExp);    
    
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
end
% -------------------------------------------------------------------------------------------------


%==================================================================================================
%% % MAKE OPTIC FLOW COMBINED VIDS------------------------------------------------------------------
%==================================================================================================

expDates = {'2018_02_07_exp_1' ...
            '2018_02_07_exp_2' ...
            '2018_02_07_exp_3' ...
            }

sids = [ 0 ...
         0 ...
         0 ...
         0 ...
         0 ...
    ];
     
 FRAME_RATE = 25;

 for iExp = 1:length(expDates)
     
     expDate = expDates{iExp}
     sid = sids(iExp)
     
     %%% Define ROIs for optic flow combined vids
     parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
     select_video_ROIs(parentDir, sid);
     writeToLog(sprintf('%s optic flow ROIs defined', expDate));
 end
 %%
 for iExp = 1:length(expDates)
     
     expDate = expDates{iExp}
     sid = sids(iExp);
     fileName = 'Behavior_Vid_ROI_Data.mat';
     
     %%% Make optic flow combined vids
     parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
     tic
     make_optic_flow_vid(sid, parentDir, FRAME_RATE, fileName)
     disp(['Creating combined vids took ', num2str(round(toc/60, 1)) ' min']);
     writeToLog(sprintf('%s combined vids created in %s min', expDate, num2str(round(toc/60, 1))));
 end

% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% ARCHIVE FILES-----------------------------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDate = '2018_01_26_exp_1';
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
trialDuration = 30;
expDate = '2018_01_17_exp_2';
sid = 3;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
annotationFileName = [expDate, '_Annotation.txt'];

tic
test = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', num2str(round(toc/60, 1)) ' min']);

% -------------------------------------------------------------------------------------------------

%% ==================================================================================================
