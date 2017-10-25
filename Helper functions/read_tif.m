function [ output ] = read_tif( tifPath )
%===================================================================================================
% Reads a .tif stack of 2P data into an array
% 
% Input: 
%       tifPath = String specifying path to .tif file
% Output: 
%       output = image data array: [Lines, Pixels, Planes, Volumes]
% ==================================================================================================

% Open .tif file
tifObj = Tiff(tifPath,'r');

% Iterate over the .tif file to find out how many images there are
nFrames = 1;
while ~tifObj.lastDirectory()
    try
        tifObj.nextDirectory();
    catch
        break;
    end
    nFrames = nFrames + 1;
end

% Extract info from image description
frameString = tifObj.getTag('ImageDescription');
nPlanes = frameStringKeyLookup(frameString, 'hFastZ.numFramesPerVolume');
nVolumes = frameStringKeyLookup(frameString, 'hFastZ.numVolumes');
nChannels = length(frameStringKeyLookup(frameString, 'scanimage.SI.hChannels.channelSave'));

% Get TIFF tag information
numLines = tifObj.getTag('ImageLength');
numPixels = tifObj.getTag('ImageWidth');
switch tifObj.getTag('SampleFormat')
    case 1
        imageDataType = 'uint16';
    case 2
        imageDataType = 'int16';
    otherwise
        assert('Unrecognized or unsupported SampleFormat tag found');
end


% Save image stack to array, accounting for multiple channels if necessary
chanData = [];
for iChannel = 1:nChannels
    tifData = zeros(numLines,numPixels,nFrames,imageDataType);  % --> [x, y, frame]
    for iFrame = 1:nFrames
        tifObj.setDirectory(iFrame);
        tifData(:,:,iFrame) = tifObj.read(); % --> [x, y, frame]
    end
    currChanData = tifData(:,:,iChannel:nChannels:end);
    chanData(:,:,:,:, iChannel) = reshape(currChanData, [numLines, numPixels, nPlanes, nVolumes]);
end

output = chanData;
tifObj.close();

end

