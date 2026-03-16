function AlignedData=AlignEyeDataToSaccades(Data,EyeTime,MarkerTime,Interval);
%This function is to align eye data according to the interval separated by
%saccades;
%Xuefei Yu;06202020
MaxCol=(diff(Interval))+1;
AlignedData=NaN*ones(numel(MarkerTime),MaxCol);

for i=1:size(MarkerTime,2)
    AlignTime=MarkerTime(:,i);
    TimeAfterAlign=repmat(EyeTime,size(MarkerTime,1),1)-AlignTime;
   for j=1:size(Data,1)
       Selection=(TimeAfterAlign(j,:)>=Interval(1)) & (TimeAfterAlign(j,:)<=Interval(2)); 
       currData=Data(j,:);
       currDataSelection=currData(Selection);
       StartBin=TimeAfterAlign(j,1)-Interval(1)+1;
       if ~isnan(StartBin)
        StartBin=max(StartBin,1);
       
          AlignedData((i-1)*size(Data,1)+j,[StartBin:StartBin+length(currDataSelection)-1])=currDataSelection;
       end
   end
    
    
    
end











end