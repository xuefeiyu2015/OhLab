function RetrivedData=RetriveDataLib(VarStr,DataLib,MeanFlag)
% Function tool to retrive variable into a matrix in the batch analysis
%MeanFlag=1: Get the first col;2: Get the second col;0:All cols
%Xuefei Yu 06192020

tmp=cellfun(@(x) x(VarStr),DataLib,'Uniform',0);
if MeanFlag
tmp=cellfun(@(x) x(:,MeanFlag)',tmp,'Uniform',0);
end

Max_Col_Whole=max(cell2mat(cellfun(@(x) size(x,2),tmp,'Uniform',0)));

RetrivedData=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),Max_Col_Whole-size(x,2))]',tmp,'Uniform',0))';


end