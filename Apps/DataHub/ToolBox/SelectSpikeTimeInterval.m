function MatrixSelected = SelectSpikeTimeInterval(Matrix_Ori,StartMat,EndMat)
%This function is for matrix data selection from start index to the end
%index for all trials; Xuefei Yu; 08112019
%Constrained to spike time; %Xuefei Yu,06082020

Matrix=Matrix_Ori(:,1:end-2);
Matrix_Start=Matrix_Ori(:,end-1);
Matrix_End=Matrix_Ori(:,end);




%   First check the dimension consistency
if (size(Matrix,1)~=size(StartMat,1))||(size(Matrix,1)~=size(EndMat,1))||(size(StartMat,1)~=size(EndMat,1))
    error('The rows of the inputs are not consistent!');
end

% Set up a matrix large enough to store all the selections 
RowNum=size(Matrix,1);
MatrixSelected_temp=ones(size(Matrix))*NaN;
maxLength=1;
if size(StartMat,2)>1
    for i=1:size(StartMat,2)
        StartInterval=StartMat(:,i);
        EndInterval=EndMat(:,i);
    
    %Check vaididty
    Check=EndInterval-StartInterval;
     NegativeTime=Check<0;
     if sum( NegativeTime)>0
        AllNormalTime=Check(~NegativeTime);
        AllNormalTime=AllNormalTime(~isnan(AllNormalTime));
       AllNormalTime_Sort=sort(AllNormalTime);
       Mid_AllNormalTime_Sort=AllNormalTime_Sort(floor(length(AllNormalTime_Sort)/2));
       EndInterval(NegativeTime)= Mid_AllNormalTime_Sort;
     end
     for j=1:RowNum
          Selection=(Matrix(j,:)>=StartInterval(j))&(Matrix(j,:)<=EndInterval(j));
           if sum(Selection)>maxLength
              maxLength=sum(Selection);
           end
           
             MatrixSelected_temp(j,1:sum(Selection))=Matrix(j,Selection);
         
     end
     
     MatrixSelected{i}=MatrixSelected_temp(:,1:maxLength);
     
    end
    
    
    
    
else

for i=1:RowNum
    Selection=(Matrix(i,:)>=StartMat(i,:))&(Matrix(i,:)<=EndMat(i,:));
    if sum(Selection)>maxLength
        maxLength=sum(Selection);
    end
    
    MatrixSelected_temp(i,1:sum(Selection))=Matrix(i,Selection);
    
end

MatrixSelected=MatrixSelected_temp(:,1:maxLength);
end

end

