function [baselineData, respData] = extract_event_volumes(eventList, filterVec, baselineDur, respDur, infoStruct, varargin)
%===========================================================================================================================
% 
% This function accepts a list of annotation events and a logical filtering vector specifying which events to extract, 
% then returns the imaging data corresponding to a baseline period and a response period for each event. If the event 
% duration is variable
%
% INPUTS:
%       eventList      = an M x 3 array specifying the onset volume, offset volume, and trial number of an event (as 
%                        columns, in that order)
%
%`      filterVec      = an M x 1 logical vector indicating which rows of eventList to use.
%
%       baselineDur    = the time in seconds before the response period to use as a baseline.
%
%       respDur        = length of the response period in seconds.
%
%       infoStruct     = structure containing the imaging data and information about it. Must have the fields
%                       [wholeSession], [nVolumes], and [volumeRate].
%
%       offsetAlign    = <OPTIONAL> a boolean value indicating whether to align the baseline and response data to the 
%                        beginning or the end of each event (default = 0).   
% OUTPUTS:
%       baselineData   = array of baseline data for each event with dimensions: [y, x, plane, volume, event] 
% 
%       respData       = same as baselineData, but for the response period instead of the baseline.
%
%=========================================================================================================================

% Parse optional arguments
p = inputParser;
addOptional(p, 'optArg', []);
addParameter(p, 'offsetAlign', 0);
parse(p, varargin{:});
offsetAlign = p.Results.offsetAlign;
if ~isempty(p.Results.optArg)
   if p.Results.optArg == 1
      offsetAlign = p.Results.optArg; % If last argument has no name/value pair, still assume it is offsetAlign if == 1
   end
end

% Set up variables
baselineDurVols = sec2vols(baselineDur, infoStruct.volumeRate);
respDurVols = sec2vols(respDur, infoStruct.volumeRate);
wholeSession = infoStruct.wholeSession;
sessionSize = size(wholeSession);
filteredList = eventList(filterVec, :);

% Calculate starting and ending volumes for each event and run validation checks
blStartVols = []; respEndVols = []; 
for iEvent = 1:size(filteredList, 1)
    
    [eventStartVol, eventEndVol, currTrial] = split_vector(filteredList(iEvent, :));
        
    % Calculate starting, ending and alignment volumes
    if offsetAlign
        alignVols(iEvent) = eventEndVol;
    else
        alignVols(iEvent) = eventStartVol;
    end
    blStartVols(iEvent) = alignVols(iEvent) - baselineDurVols - 1; % Subtract one more so the baseline ends just before stimEndVol
    respEndVols(iEvent) = alignVols(iEvent) + respDurVols;
        
end

alignVols(blStartVols == 0) = alignVols(blStartVols == 0) + 1;
respEndVols(blStartVols == 0) = respEndVols(blStartVols == 0) + 1;
blStartVols(blStartVols == 0) = 1;


% Loop through and pull out data
baselineData = zeros([sessionSize(1:3), baselineDurVols + 1, size(filteredList, 1)]);
respData = zeros([sessionSize(1:3), respDurVols + 1, size(filteredList, 1)]);
for iEvent = 1:size(filteredList, 1)
    currTrial = filteredList(iEvent, 3);
    baselineData(:,:,:,:,iEvent) = wholeSession(:,:,:, blStartVols(iEvent):(alignVols(iEvent)-1), currTrial);    % --> [y, x, plane, volume, event]
    respData(:,:,:,:,iEvent) = wholeSession(:,:,:, alignVols(iEvent):respEndVols(iEvent), currTrial);            % --> [y, x, plane, volume, event]
end

end