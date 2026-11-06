function [ColorHigh, ColorLow] = ValueColors(Context)
%VALUECOLORS Lab colour convention for high / low value.
%   Object value   : high = red,    low = blue
%   Location value : high = purple, low = green
%
% Input:
%   Context : 'Object' / 'Obj' / 'Location' / 'Loc', or a task code
%             (ONEDR_TASK and OBJ_LOCATIONVALUE_TASK -> location; any
%             other task code -> object)
% Outputs:
%   ColorHigh, ColorLow : 1x3 RGB
%
% Xuefei Yu

ECode;
if isnumeric(Context)
    IsLocation = ismember(Context, [ONEDR_TASK, OBJ_LOCATIONVALUE_TASK]);
else
    IsLocation = any(strcmpi(Context, {'Location', 'Loc'}));
end

if IsLocation
    ColorHigh = [0.49 0.18 0.56];   % purple
    ColorLow  = [0.13 0.69 0.47];   % green
else
    ColorHigh = [0.85 0.15 0.15];   % red
    ColorLow  = [0.15 0.30 0.85];   % blue
end
end
