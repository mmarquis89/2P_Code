function nFrames = count_frames(vidFile)
%===================================================================================================
% Loads a video file and returns the number of frames in it.
% 
% Inputs:
%       vidFile = the path (incl. file name) to the file you want to process
%===================================================================================================
    
    myVid = VideoReader(vidFile);
    nFrames = 0;
    while hasFrame(myVid)
        readFrame(myVid);
        nFrames = nFrames + 1;
    end
end