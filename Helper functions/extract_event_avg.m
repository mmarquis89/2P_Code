function [baselineData, eventData] = extract_event_avg(eventList, filterVec, baselineDur, infoStruct, varargin)
%=========================================================================================================================
% 
% This function accepts a list of annotation events and a logical filtering vector specifying which events to extract, 
% then returns the imaging data corresponding to a baseline period and the event itself, averaged across volumes. The 
% baseline period can be of any length and can be either before or after the event.
%
% INPUTS:
%       eventList    = an M x 3 array specifying the onset volume, offset volume, and trial number of an event (as 
%                      columns, in that order)
%
%`      filterVec    = an M x 1 logical vector indicating which rows of eventList to use.
%
%       baselineDur  = the time in seconds of the baseline period.
%
%       infoStruct   = structure containing the imaging data and information about it. Must have the fields
%                     [wholeSession], [nVolumes], and [volumeRate].
%
%       offsetAlign  = <OPTIONAL> a boolean value indicating whether to collect the baseline data from the  
%                      beginning or the end of each event (default = 0).   
%
% OUTPUTS:
%       baselineData = array of averaged baseline data for each event with dimensions: [y, x, plane, event] 
% 
%       eventData    = same as baselineData, but for the volumes of the event itself.
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
baselineDurVols = floor(baselineDur * infoStruct.volumeRate);
wholeSession = infoStruct.wholeSession;
sessionSize = size(wholeSession);
filteredList = eventList(filterVec, :);

% 
baselineVols = []; eventVols = [];
for iEvent = 1:size(filteredList, 1)
    
    % Get list of current event volumes
    [eventStartVol, eventEndVol, ~] = split_vector(filteredList(iEvent, :));    
    eventVols{iEvent} = eventStartVol:eventEndVol;
    
    % Calculate baseline volumes
    if offsetAlign
        baselineVols{iEvent} = (eventEndVol + 1):(eventEndVol + baselineDurVols + 1);
    else
        baselineVols{iEvent} = (eventStartVol - baselineDurVols - 1):(eventStartVol - 1);
    end
    
    % Be sure not to try and get volumes before the start or after the end of the trial
    if baselineVols{iEvent}(1) < 1
        baselineVols{iEvent} = 1:eventStartVol - 1;
    end
    if baselineVols{iEvent}(end) > infoStruct.nVolumes
        baselineVols{iEvent} = (eventEndVol + 1):infoStruct.nVolumes;
    end
end

% Loop through and pull out data
baselineData = zeros([sessionSize(1:3), size(filteredList, 1)]); eventData = baselineData;
for iEvent = 1:size(filteredList, 1)
    currTrial = filteredList(iEvent, 3);
    eventData(:,:,:,iEvent) = squeeze(mean(wholeSession(:,:,:, eventVols{iEvent}, currTrial), 4));          % --> [y, x, plane, event]
    baselineData(:,:,:,iEvent) = squeeze(mean(wholeSession(:,:,:, baselineVols{iEvent}, currTrial), 4));    % --> [y, x, plane, event]
end

end