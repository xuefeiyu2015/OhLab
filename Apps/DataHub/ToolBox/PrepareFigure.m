function fighandle=PrepareFigure(Num,Color,Position,varargin)
%General procedural to pre-prepare the figure
%Xuefei Yu


    if ishandle(Num)
    close(Num);end


fighandle=figure(Num);
set(Num,'color',Color);
set(Num,'Position', Position);

if ~isempty(varargin) && contains(varargin{1},'Name')
    set(Num,'Name',varargin{2});
end


end