function mat=fillinNaNs(cell_ori)
%Function to convert cell to matrix by filling in NaN at the end
%Xuefei Yu; 06192020
maxnumel=max(cellfun(@(x) size(x,2),cell_ori));
maxrow=max(cellfun(@(x) size(x,1),cell_ori));
mat=cell2mat(cellfun(@(x) [x, NaN*ones(max(maxrow-size(x,1),size(x,1)),maxnumel-size(x,2))],cell_ori,'uniform',0));


end
