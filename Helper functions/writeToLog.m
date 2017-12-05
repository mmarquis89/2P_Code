function writeToLog(msg)

% Appends lines to a log file for overnight data processing
myFile = fopen('D:\Dropbox (HMS)\2P Data\Imaging Data\dataProcessingLog.txt', 'a');
fprintf(myFile, [datestr(datetime), ':  ', msg '\r\n']); 
fclose(myFile);

end