function [ output ] = readTif( tifPath )
% Input: String with path to .tif file
% Output: [Lines, Pixels, Planes, Volumes]

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

% Save image stack to array
output = zeros(numLines,numPixels,nFrames,imageDataType);
for iframe = 1:nFrames
    tifObj.setDirectory(iframe);
    output(:,:,iframe) = tifObj.read();
end

output = reshape(output, [numLines,numPixels,nPlanes,nVolumes]);

tifObj.close();

end

