function T = updateTable(T,varargin);
%
% T = updateTable(T, var1, var2, ...)
%
% Adds a new row to the table T consisting of the input variables.
% 
% The total number of vars (var1, var2, ... varN) must be equal to
% the number of columns in T.
%
% example of usage:
% T = initTable('conditions','values');
% T = updateTable(T,1,100);
% T = updateTable(T,2,101);
%
% last modified 2016-dec-06
% dbtm

cols = T.Properties.VariableNames;
if length(cols)~=length(varargin)
    error('Number of vars must equal the number of columns in T.');
    help updateTable
end

argString = '(';
for n=1:length(cols)
    
    eval([cols{n} ' = [' num2str(varargin{n}) '];']);
    argString = [argString cols{n}];
    if n<length(cols)
        argString = [argString ','];
    end
end
argString = [argString ')'];

newRow = eval(['table' argString ';']);
T = vertcat(T,newRow);
