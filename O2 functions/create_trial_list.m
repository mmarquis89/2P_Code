function create_trial_list(parentDir, sid)

% Create a comma-separated text file listing all the folders to be processed 
% for a specific session and another one listing the tid of each dir

disp(parentDir)

dirContents = dir(fullfile(parentDir, '*sid*tid*'));
dirFile = fopen(fullfile(parentDir, 'dirList.txt'), 'w');
disp(num2str(numel(dirContents)))
for i = 1:numel(dirContents)
    disp(fullfile(parentDir, dirContents(i).name))
    if isdir(fullfile(parentDir, dirContents(i).name))
        currSid = str2double(regexp(dirContents(i).name, '(?<=sid_).*(?=_tid)', 'match'));
        if currSid == sid
            tid = regexp(dirContents(i).name, '(?<=tid_).*(?=_)', 'match');
            writeStr = [fullfile(parentDir, [dirContents(i).name]), ',', tid{:}, '\n'];
            fprintf(dirFile, writeStr);
        end
    end
end
fclose(dirFile);