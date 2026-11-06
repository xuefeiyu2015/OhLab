function Area = ComputeBatchRowArea(FileName, ChannelNumber, Ref)
% ComputeBatchRowArea  Assign each batch row (session + channel range) to a
% brain area by matching it against reference batch lists. Matching uses the
% session name AND channel range, because simultaneous FEF+SC sessions
% (e.g. Adams042624FRACSAC1) appear in both lists with different channels.
%
% Inputs
%   FileName      — cell{nRow} session names (Data.FileName)
%   ChannelNumber — nRow-by-2 [startChan endChan] (Data.ChannelNumber)
%   Ref           — struct array with fields .Area (char), .FileName (cell),
%                   .ChannelNumber (nRef-by-2), one per reference batch file
% Output
%   Area — cell{nRow} area label per row ('' if no reference row matches)
%
% Xuefei Yu, 2026-11-06

nRow = numel(FileName);
Area = repmat({''}, nRow, 1);
for i = 1:nRow
    for r = 1:numel(Ref)
        hit = strcmp(Ref(r).FileName(:), FileName{i}) & ...
              all(Ref(r).ChannelNumber == ChannelNumber(i,:), 2);
        if any(hit)
            Area{i} = Ref(r).Area;
            break;
        end
    end
end
end
