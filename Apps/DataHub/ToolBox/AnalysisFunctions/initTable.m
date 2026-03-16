function T = initTable(varargin);
%
% T = initTable(columName1, columName2, ...);
%
% Creates a blank table with zero rows and N columns, where N is
% the number of input arguments. 
% 
% Input arguments must be character strings specifying column
% names.
%
% example of usage:
% T = initTable('conditions','values');
% T = updateTable(T,1,100);
% T = updateTable(T,2,101);
%
% last modified 2016-dec-06
%
argString = '(';
for n=1:length(varargin)
    eval([varargin{n} ' = [];']);
    argString = [argString varargin{n}];
    if n<length(varargin)
        argString = [argString ','];
    end
end
argString = [argString ')'];

T = eval(['table' argString ';']);
