function OrganizedData = ReOrganizeData(Data)
%Reorganize the data according to the trials;
%Each line is one trial
%YXF 10142019

if isempty(Data)
    OrganizedData = Data;
    return;
end
StackNum=max((size(Data,2)-1),1);
LineNum=unique(Data(:,1));





RepNum = countmember(unique(Data(:,1)),Data(:,1));
TrialSeq=unique(Data(:,1));

if (length(unique(RepNum))==1 && unique(RepNum)==1)
    %If the orginial data set is already organized as required,we can save the
    %time
     OrganizedData=squeeze(Data);
     return
else
    RepNum=max(RepNum)+1;
    LineNum=length(LineNum);

%Initiate the organized data matrix
    OrganizedData=NaN*zeros(LineNum,RepNum,StackNum);
    for i=1:LineNum
        for j=1:StackNum
        select=Data(Data(:,1)==TrialSeq(i),j+1);
        OrganizedData(i,1:(numel(select)+1),j)=[TrialSeq(i),reshape(select,1,numel(select))];
  
        end
 
    end
 
    OrganizedData=squeeze(OrganizedData);
end

end

