function range = calc_range(dataArr, scalar)
%======================================================================================================
% Calculates the maximum absolute value of an array and returns [min, max] values for the array that
% are symmetrically centered at zero. Useful for example in calculating ranges for plots using a 
% colormap centered at zero. Can optionally multiply output by a scalar value.
%
% INPUTS:
%       dataArr = numeric array of any dimension containing the data to be encompassed by the range.
%       scalar  = <optional> scalar value to multiply the output range by. Can pass [] to omit scaling.
%
% OUTPUTS:
%       range   = numeric vector in the form of [min, max] containing the range values (e.g., [-3, 3]) 
%=======================================================================================================

if isempty(scalar)
    scalar = 1;    
end

minMax = [];
minMax(1) = min(dataArr(:));
minMax(2) = max(dataArr(:));
range(2) = max([abs(min(minMax)), max(minMax)]);
range(1) = -range(2);
range = range * scalar;

end