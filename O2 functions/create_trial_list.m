function create_trial_list(parentDir)

% Create a comma-separated text file listing all the folders to be processed, 
% and another one listing the sid and tid of each dir

disp(parentDir)

dirContents = dir(fullfile(parentDir, '*sid*tid*'));
dirFile = fopen(fullfile(parentDir, 'dirList.txt'), 'w');
disp(num2str(numel(dirContents)))
for i = 1:numel(dirContents)
    disp(fullfile(parentDir, dirContents(i).name))
    if isdir(fullfile(parentDir, dirContents(i).name))
        sid = regexp(dirContents(i).name, '(?<=sid_).*(?=_tid)', 'match');
        tid = regexp(dirContents(i).name, '(?<=tid_).*(?=_)', 'match');
        name = dirContents(i).name;
        test = whos('name');
        disp(test.class)
        writeStr = [fullfile(parentDir, [dirContents(i).name]), ',', sid{:}, ',', tid{:}, '\n'];
        fprintf(dirFile, writeStr);
    end
end
fclose(dirFile);