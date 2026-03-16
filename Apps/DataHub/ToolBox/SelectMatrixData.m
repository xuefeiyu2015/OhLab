function MatrixSelected = SelectMatrixData(Matrix,StartMat,EndMat)
%This function is for matrix data selection from start index to the end
%index for all trials; Xuefei Yu; 08112019
%   First check the dimension consistency
if (size(Matrix,1)~=size(StartMat,1))||(size(Matrix,1)~=size(EndMat,1))||(size(StartMat,1)~=size(EndMat,1))
    error('The rows of the inputs are not consistent!');
end
% Set up a matrix large enough to store all the selections 
Col_Max=max(EndMat-StartMat+1);
RowNum=size(Matrix,1);

MatrixSelected=ones(RowNum,Col_Max)*NaN;


for i=1:RowNum
    SelectLength=EndMat(i)-StartMat(i)+1;
    CurrentMat=Matrix(i,:);
    if StartMat(i)<0
        MatrixSelected(i,1:-StartMat(i)+1)=NaN;
        MatrixSelected(i,-StartMat(i)+1+1:EndMat(i)-StartMat(i)+1)=CurrentMat(1:EndMat(i));
        
    else
    
    MatrixSelected(i,1:SelectLength)=CurrentMat(StartMat(i):EndMat(i));
    end
    
end

end

