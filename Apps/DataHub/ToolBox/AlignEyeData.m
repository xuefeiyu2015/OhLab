function AlignedData=AlignEyeData(Data,EyeTime,MarkerTime,Interval);
%This function is to align eye data according to the interval separated by
%markers ;
%Xuefei Yu;05122020
MaxCol=ceil(Interval(2)-Interval(1))+1;
AlignedData=NaN*ones(size(Data,1),MaxCol);


for i=1:size(Data,1)
    TimeAfterAlign=EyeTime-MarkerTime(i);
    currInterVal=Interval;
    Selection=(TimeAfterAlign>=currInterVal(1)) & (TimeAfterAlign<=currInterVal(2)); 
    currData=Data(i,:);
    currDataSelection=currData(Selection);
    StartBin=TimeAfterAlign(1)-currInterVal(1)+1;
    if ~isnan(StartBin)
        StartBin=max(StartBin,1);
    
        AlignedData(i,[StartBin:StartBin+length(currDataSelection)-1])=currDataSelection;
    
        
    end
    
end







end