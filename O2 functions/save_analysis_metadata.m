function initial_analysis_processing(parentDir, sessionDataFile)

    
    % Create and save a structure of analyis metadata + some hardcoded parameters
    [analysisMetadata, ~] = load_imaging_data(parentDir, sessionDataFile);
    analysisMetadata.ROIdata = [];
    analysisMetadata.MAX_INTENSITY = 800; % To control brightness of ref image plots
    analysisMetadata.FRAME_RATE = 25; % This is the frame rate of the behavior video, not the GCaMP imaging
    save(fullfile(parentDir, 'analysisMetadata.mat'), 'analysisMetadata')
    
    % Process and save annotation types
    skipTrials = [];
    if exist(fullfile(parentDir, 'skipTrials.mat'), 'file')
        load(fullfile(parentDir, 'skipTrials.mat')) % should contain one variable called 'skipTrials'
    elseif exist(fullfile(parentDir, 'skipTrials.txt'), 'file')
        myFile = fopen(fullfile(parentDir, 'skipTrials.txt', 'r'));
        currLine = fgetl(myFile);
        while ischar(currLine)
            skipTrials(end + 1) = str2double(currLine);
            currLine = fgetl(myFile);
        end
        fclose(myFile);      
    end
    
    
    
    
    
end