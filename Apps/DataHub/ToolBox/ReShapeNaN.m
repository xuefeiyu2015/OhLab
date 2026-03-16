function reshapedMat=ReShapeNaN(mat,ord)
%This function is to reshape matrix to put all nan in the front or in the
%back
%,the second varag; 1:put all nan in the front;-1:put all nan in the back
%Xuefei Yu 03032020
reshapedMat=NaN*ones(size(mat));
for i=1:size(mat,1)
    orimat=mat(i,:);
    orimat=orimat(~isnan(orimat));
    if ord==1
    reshapedMat(i,length(reshapedMat(i,:))-length(orimat)+1:end)=orimat;
    elseif ord==-1
    reshapedMat(i,1:length(orimat))=orimat;
    end
end




end