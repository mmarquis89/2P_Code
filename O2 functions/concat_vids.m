function concat_vids(vidDir, fileStr, varargin)
%===================================================================================================
% CONCATENATE A SERIES OF VIDEOS INTO A SINGLE FILE
% Concatenates a series of videos (e.g. for each trial of a 2P experiment). The new 
% video will be saved in the same folder as the source videos, along with a count of the total
% number of frames in the concatenated video.
% 
% NOTE: the file names will be sorted in alphabetical order to determine sequence of concatenation.
%
% INPUTS:
% 
%       vidDir = the directory containing the videos you want to combine 
%                   e.g. 'U:\2P Behavior Video\2017_07_30\_Movies'
% 
%       fileStr = a filtering string to identify the files to be processed, e.g. '*sid_0_tid*.mp4'. 
%                 Be careful to ensure that only the desired files will meet this specification.
% 
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%
%       'OutputDir'  = (default: vidDir) the directory to save the output video in.
%
%       'FrameRate'  = (default: 25) the frame rate that the video was acquired at in FPS
%
%       'OutputFile' = (default = fileStr + '_AllTrials') the desired name of the output file (minus
%                      the file extension). Default removes all wildcard characters in the fileStr 
%                      and appends '_AllTrials'
%
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'FrameRate', 25);
addParameter(p, 'OutputFile', [regexprep(regexprep(regexprep(fileStr, '*', ''), '.mp4', ''), '.avi', ''), '_AllTrials']);
parse(p, varargin{:});
FRAME_RATE = p.Results.FrameRate;
outputFileName = p.Results.OutputFile;

% Identify files to be processed
vidFiles = dir(fullfile(vidDir, fileStr));
vidNames = sort({vidFiles.name});

% Create vidWriter
myVidWriter = VideoWriter(fullfile(vidDir, outputFileName), 'Uncompressed AVI');
myVidWriter.FrameRate = FRAME_RATE;
open(myVidWriter)
frameCount = 0;

disp('Concatenating videos...')
for iTrial = 1:length(vidNames)
       
    disp(vidNames{iTrial})
            
        % Load movie for the current trial
        myMovie = {};
        myVid = VideoReader(fullfile(vidDir, vidNames{iTrial}));
        while hasFrame(myVid)
            currFrame = readFrame(myVid);
            myMovie(end+1) = {uint8(currFrame)};
        end
        
        % Add frames to movie
        for iFrame = 1:length(myMovie)
            writeVideo(myVidWriter, myMovie{iFrame});
            frameCount = frameCount + 1;
        end
end%iTrial

close(myVidWriter)
clear('myMovie')

% Save frame count log
save(fullfile(vidDir, [outputFileName, '_frameCountLog.mat']), 'frameCount')

end%function