function extract_ROI_data(parentDir, sessionDataFile, ROIfile)

% Load ROI info
load(fullfile(parentDir, ROIfile)); % Contains variable 'ROIdata'
nROIs = numel(ROIdata);
disp('ROIs loaded')

% Load imaging data
[analysisMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile);

% Extract and save xy-averaged data from each ROI
ROIDataAvg = [];
for iROI = 1:nROIs
    
    disp(['Extracting data for ROI #', num2str(iROI), ' of ', num2str(nROIs), '...'])
    currMask = ROIdata(iROI).mask;
    currPlane = ROIdata(iROI).plane;
    currPlaneData = squeeze(wholeSession(:,:,currPlane,:,:));                                   % --> [y, x, volume, trial]
    currPlaneData(~currMask(:,:,ones(1, nVolumes), ones(1, nTrials))) = nan;                    % --> [y, x, volume, trial]
    currDataLin = reshape(currPlaneData, size(currPlaneData, 1)*size(currPlaneData, 2), ...
        analysisMetadata.nVolumes, analysisMetadata.nTrials);                                   % --> [pixel, volume, trial, ROI]
    ROIDataAvg(:,:,iROI) = squeeze(mean(currDataLin, 1, 'omitnan'));                            % --> [volume, trial, ROI]
end
save(fullfile(parentDir, 'ROI_Data_Avg.mat'))

end