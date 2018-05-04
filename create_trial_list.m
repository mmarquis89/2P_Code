function create_trial_list(expDate)

% Create a comma-separated text file listing all the folders to be processed, 
% and another one listing the sid and tid of each dir
parentDir = fullfile('/n/scratch2', expDate, 'BehaviorVideo');
dirContents = dir(fullfile(parentDir, '*sid*tid*'));
dirFile = fopen(fullfile(parentDir, 'dirList.txt'), 'w');
for i = 1:numel(dirContents)
    if isdir(fullfile(parentDir, dirContents(i).name))
        sid= regexp(dirContents(i).name, '(?<=sid_).*(?=_)', 'match');
        tid = regexp(dirContents(i).name, '(?<=tid_).*(?=_)', 'match');
        fprintf(dirFile, [fullfile(parentDir, dirContents(i).name), ',', sid, ',', tid, '\n']);
    end
end