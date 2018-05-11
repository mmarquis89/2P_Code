function make_vid(inputDir, sid, tid, varargin)
%==============================================================================================================================
% CREATE MOVIES FROM .TIF FILES
% Creates a .mp4 movie from a directory of .tif files captured by the fly behavior camera.
%
% Inputs:
%
%       inputDir = the directory containing the video frames you want to create a movie from.
%                   e.g. 'U:\2P Behavior Video\2017_07_30'
%
%       sid       = the session ID of the video you want to process.
%
%       tid       = the trial ID of the video you want to process.
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       OutputDir        = (default: dirPath) directory to save the output files in
%
%       FrameRate        = (default: 25) the frame rate that the behavior camera was acquiring at.
%
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutputDir', inputDir);
addParameter(p, 'FrameRate', 25);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;
frameRate = p.Results.FrameRate;

frameFiles = dir(fullfile(inputDir, 'fc2*.tif'));

% Pad the trial number with leading zeros if necessary to ensure correct filename sorting
if tid < 10
    padStr = '00';
elseif tid < 100
    padStr = '0';
else
    padStr = '';
end

% Initialize variables
trialStr = ['sid_', num2str(sid), '_tid_', padStr, num2str(tid)];
% Create save directory if it doesn't already exist
if ~isdir(outputDir)
    mkdir(outputDir);
end

% Make sure this video doesn't already exist
if exist(fullfile(outputDir, [trialStr, '.mp4']), 'file') == 0
    
    % Create video writer object
    outputVid = VideoWriter(fullfile(outputDir, trialStr), 'Motion JPEG AVI');
    outputVid.FrameRate = frameRate;
    open(outputVid)
    
    % Make sure there's at least one image file in this trial's directory
    if ~isempty(frameFiles)
        currFrames = sort({frameFiles.name}');
        
        % Write each .tif file to video
        for iFrame = 1:length(currFrames)
            
            % Read image
            currImg = imread(fullfile(inputDir, currFrames{iFrame}));
            
            % Write frame to video
            writeVideo(outputVid, currImg);
        end
    end%if
    close(outputVid)
else
    disp(['Video already exists...skipping ', trialStr]);
end%if

end%function