%==================================================================================================
%%% EXTRACT SAMPLE FRAMES FROM ANATOMY STACKS -----------------------------------------------------
%% ==================================================================================================

expDates = {...
    '2018_04_16'
            }

targetPlanes = [ 200 ...
    ];

fileStr = '*Stack_*.tif';

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    targetPlane = targetPlanes(iExp);
    
    dirPath = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    disp(['Extracting sample frames from anatomy stacks...']);
    extract_sample_frames(dirPath, fileStr, targetPlane);
    writeToLog(sprintf('%s sample frames extracted in %s min', expDate, num2str(round(toc/60, 1))));
    disp(['Extracting sample frames took ', num2str(round(toc/60, 1)) ' min']);
end
clear expDates targetPlanes fileStr dirPath
% -------------------------------------------------------------------------------------------------

%% ===================================================================================================
%%% SELECT ROI FOR OPTIC FLOW CALCULATION
%% ===================================================================================================

expDates = {...
    '2018_04_16'
            }

sids = [ ...
    0 ...
    ];
     

 FRAME_RATE = 25;

 for iExp = 1:length(expDates)
     
     expDate = expDates{iExp}
     sid = sids(iExp)
     
     %%% Define ROIs for optic flow combined vids
     parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
     select_video_ROIs(parentDir, sid);
     writeToLog(sprintf('%s optic flow ROIs defined', expDate));
 end


%% ==================================================================================================
%%% CREATE AVERAGE FLUORESCENCE VIDEOS ------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDates = {...
    '2018_04_16'
            }

sids = [ ...
    0 ...
    ];

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    sid = sids(iExp);
    
    tic
    parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    disp('Creating average fluorescence videos...');
    create_average_fluorescence_vid(parentDir, sid)
    writeToLog(sprintf('%s average fluorescence vids made in %s min', expDate, num2str(round(toc/60, 1))));
    disp(['Creating average fluorescence videos took ', num2str(round(toc/60, 1)) ' min']);
end
clear expDates sids parentDir
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% CREATE ANATOMY STACKS --------------------------------------------------------------------------
%==================================================================================================

expDates = {...
    '2018_04_16'
            }

fileStr = '*Stack*.tif';

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    
    parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    create_anatomy_stack(parentDir, fileStr, '_Stack');
    disp(['Creating the anatomy stacks took ', num2str(round(toc/60, 1)) ' min']);
    writeToLog(sprintf('%s anatomy stack created in %s min', expDate, num2str(round(toc/60, 1))));
end
clear expDates fileStr parentDir
% -------------------------------------------------------------------------------------------------

%===================================================================================================
%%% PRE-REGISTRATION PROCESSING -------------------------------------------------------------------
%===================================================================================================
% 
% <3 you Michael!
% -------------------------------------------------------------------------------------------------

expDates = {...
    '2018_04_16'
            }

sids = [ ...
    0 ...
    ];

disp('Starting pre-processing...')
for iExp = 1:length(expDates)
    
    expDate = expDates{iExp};
    sid = sids(iExp);
    
    parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    preReg_routine_MM(parentDir, sid, expDate) 
    writeToLog(sprintf('%s_%s pre-processing completed in %s min', expDate, num2str(round(toc/60, 1))));
end
clear expDates sids parentDir
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% NoRMCorre REGISTRATION --------------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDates = {...
    '2018_04_16'
            }

sids = [ ...
    0 ...
    ];

% fileNames = repmat({['sid_', num2str(sid), '_sessionFile.mat']}, 1, numel(expDates));
% fileNames = { 'sid_1_Chan_1_sessionFile.mat', ...
%             }
disp('Starting NoRMCorre registration...')
for iExp = 1:length(expDates)
    
    expDate = expDates{iExp};
    sid = sids(iExp);
    fileName = ['sid_', num2str(sid), '_sessionFile.mat'];%fileNames{iExp};
    
    parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
    
    tic; normcorre_registration(parentDir, fileName);
    writeToLog(sprintf('%s_%s NoRMCorre registration completed in %s min', expDate, fileName, num2str(round(toc/60, 1))));
    
end
clear expDates sids fileNames fileName parentDir
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%%% CREATE BEHAVIOR VIDEOS ------------------------------------------------------------------------
%==================================================================================================
expDates = {...
    '2018_04_16'
            }

sids = [ ...
    0 ...
    ];

FRAME_RATE = 25;

disp('Creating behavior vids...')
for iExp = 1:length(expDates)
    
    expDate = expDates{iExp};
    sid = sids(iExp);    
    
    %%% Make vids
    parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
    tic
    make_vids_2P(sid, parentDir, FRAME_RATE)
    disp(['Creating behavior vids took ', num2str(round(toc/60, 1)) ' min']);
    writeToLog(sprintf('%s behavior vids created in %s min', expDate, num2str(round(toc/60, 1))));
    
    %%% Concatenate vids
    tic
    parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
    concat_vids_2P(sid, parentDir, FRAME_RATE)
    disp(['Concatenating behavior vids took ', num2str(round(toc/60, 1)) ' min']);
    writeToLog(sprintf('%s behavior vids concatenated in %s min', expDate, num2str(round(toc/60, 1))));
end
clear expDates sids parentDir
% -------------------------------------------------------------------------------------------------


%==================================================================================================
%%% MAKE OPTIC FLOW COMBINED VIDS------------------------------------------------------------------
%==================================================================================================
 
expDates = {...
    '2018_04_16'
            }

sids = [ ...
    0 ...
    ];

fileNames = { ...
    'Behavior_Vid_ROI_Data.mat' ...
    'Behavior_Vid_ROI_Data.mat' ...
};

disp('Creating optic flow combined vids...')
for iExp = 1:length(expDates)
     
     expDate = expDates{iExp}
     sid = sids(iExp);
     fileName = fileNames{iExp};
     
     %%% Make optic flow combined vids
     parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
     tic
     make_optic_flow_vid(sid, parentDir, FRAME_RATE, fileName)
     disp(['Creating combined vids took ', num2str(round(toc/60, 1)) ' min']);
     writeToLog(sprintf('%s combined vids created in %s min', expDate, num2str(round(toc/60, 1))));
 end
 clear expDates sids parentDir fileName
% -------------------------------------------------------------------------------------------------

%==================================================================================================
%% ARCHIVE FILES-----------------------------------------------------------------------------------
%==================================================================================================

% -------------------------------------------------------------------------------------------------

expDate = '2018_02_09_exp_1';
parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
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
parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
archiveName = ['sid_', num2str(sid), '_RawFrames'];
system7zip(parentDir, archiveName, '7z', ['*sid_', num2str(sid), '_t*'], 1);

clear parentDir archiveName filterString
% -------------------------------------------------------------------------------------------------
%==================================================================================================
%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
%==================================================================================================

FRAME_RATE = 25;
trialDuration = 20;
expDate = '2018_04_14_exp_2';
sid = 0;

parentDir = ['B:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
annotationFileName = [expDate, '_Annotation.txt'];
% annotationFileName = [expDate, '_sid_', num2str(sid), '_Annotation.txt'];

tic
process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', num2str(round(toc/60, 1)) ' min']);

clear parentDir saveDir annotationFileName
% -------------------------------------------------------------------------------------------------

%% ==================================================================================================
