function planeMovie=generatePlaneMovie_2P_akm(planeNum,path2file)
% Use output with implaymovie:
%    planeMovie=generatePlaneMovie_2P_akm(planeNum,path2file);
%    implay(planeMovie);
% In the implay GUI, open Tools > Colormap > "specify range" & enter min/max
% values for video playback, e.g 50-255.

cd(path2file)
in=load('sessionOutfile_Reg1.mat');
wholeSession=in.regProduct;
wholeSession = testSession; % [y, x ,z, volume, trialNum]

totalFrames=size(wholeSession,4)*size(wholeSession,5); % number of volumes per trial * number of trials
planeMovie=uint16(zeros( [size(wholeSession,1),size(wholeSession,1),totalFrames] )); % Initialize a new 3d array, [x y frame]

frameCount=1;
for t=1:size(wholeSession,4) % For each volume (timepoint) in each trial...
    for r=1:size(wholeSession,5)
        planeMovie(:,:,frameCount)=squeeze(wholeSession(:,:,planeNum,t,r));
        frameCount=frameCount+1;
    end
end

planeMovie2 = planeMovie(:,:,1:6:end);
planeMovie3 = planeMovie(:,:,3:6:end);

end

