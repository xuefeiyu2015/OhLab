function cmd = unpackHeader(dat);
%
% cmd = unpackHeader(dat);
%
% Utility script for unpacking the contents of dat.h.vv into the current
% memory space.  Replaces the set_[experiment]_global usages from earlier
% versions of hanuman.
% 
% INPUT: hanuman data structure.  2012-dec-31
% dbtm

names = fieldnames(dat.h.vv);

cmd = [];
for n = 1:length(names)
    val = getfield(dat.h.vv,names{n});
    if length(val)>1
        s = [names{n} ' = [' num2str(val) ']' ];
    else
        s = [names{n} ' = ' num2str(val)];
    end
    
    cmd = [cmd s ';'];
end