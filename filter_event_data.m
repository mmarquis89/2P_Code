function filtOutput = filter_event_data(eventList, filterEventData, filterEventWindows, filterDirections)
%==============================================================================================================
% 
% This function takes a list of annotation events (e.g. odor delivery, ball stopping, locomotion) and 
% determines whether they meet a list of criteria based on whether other types of events ("filter" events)
% occurred just before, during, or after the target event. The number of filter event types, size of the 
% evaluation time windows, and direction of the filtering (i.e. whether the filter event is required to be 
% present or absent to match the filter) are all adjustable.
%
% INPUTS:
%       eventList           = an M x 3 array specifying the onset volume, offset volume, and trial number of 
%                             an event (as columns, in that order)
%
%       filterEventData     = an A x B x N array containing the annotation data (in volumes) for each type of 
%                             event that you want to use to filter the target events in eventList. The 
%                             dimensions of the array should be: [trial, volume, eventType]
% 
%       filterEventWindows  = an N x 2 array specifying the [pre-stim, post-stim] time windows in seconds 
%                             around the target event onset and offset to use for filtering. Should have one 
%                             row for each type of filter event.
% 
%       filterDirections    = an N x 3 array of numeric values indicating whether each filter type should be 
%                             inclusive or exclusive during the pre-stim, during stim, and post-stim filter 
%                             windows. A value of 1 means inclusive, -1 means exclusive, and 0 means that you
%                             don't care whether the filter event happened during that period. For example, 
%                             passing [-1 -1 -1 ; -1 1 0] would match the target events for which:
%                                - The first type of filter event didn't occur at any point in the time window
%                                - The second type of filter event occurred during the target event, but not 
%                                  in the window before the event (and the post-event window is irrelevant).
%
% OUTPUTS:
%       filtOutput          = a logical vector indicating whether each event in the eventList passed the 
%                             filtering process for all the filter event types.
% 
%==============================================================================================================

% Calculate constants
nVolumes = size(filterEventData, 2);
nEvents = size(eventList, 1);
nFilters = size(filterEventData, 3);

% Fill in a logical array for each event and filter
filtArr = zeros(nEvents, nFilters);
for iFilt = 1:nFilters
    filtWin = filterEventWindows(iFilt, :);
    for iEvent = 1:nEvents
        [eventStart, eventEnd, trial] = split_vector(eventList(iEvent, :));
        
        % Calculate the starting and ending volumes of the filtering window
        startVol = eventStart - filtWin(1);
        endVol = eventEnd + filtWin(2);
        if startVol < 1
            startVol = 1;
        end
        if endVol > nVolumes
            endVol = nVolumes;
        end
        
        % Determine whether any filter events occurred within the filter window
        if sum(filterEventData(trial, startVol:(eventStart - 1), iFilt))
            preFilt = 1;
        else
            preFilt = -1;
        end
        if sum(filterEventData(trial, eventStart:eventEnd, iFilt))
            stimFilt = 1;
        else
            stimFilt = -1;
        end
        if sum(filterEventData(trial, (eventEnd + 1):endVol, iFilt))
            postFilt = 1;
        else
            postFilt = -1;
        end
        
        % Check whether the current event matches the filter
        filterVec = [preFilt, stimFilt, postFilt];
        filterVec(filterDirections(iFilt, :) == 0) = 0; % Force matching for windows you don't care about
        if sum(filterVec == filterDirections(iFilt, :)) == 3
            filtArr(iEvent, iFilt) = 1;
        end
    end
end
filtOutput = logical(sum(filtArr, 2) == nFilters);
end