function EyeInterval=SelectEyeInterval(EyeData,Start,End)
%This function is to select our the eye interval according to the Start
%matrix and the end matrix
%Xuefei Yu 03220202

for i=1: size(Start,2)
    StartInterval=round(Start(:,i));
    EndInterval=round(End(:,i));
    
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
    
    
    %Judge whether the eye data is empty
    EmptyCheck=~cellfun(@(x) sum(isnan(x))==numel(x),EyeData);
    TrialIndex_tmp=(1:length(StartInterval))';
    SelectionMat=num2cell((1:length(TrialIndex_tmp(EmptyCheck)))'); 
    
    EyeInterval_temp=cell(length(StartInterval),1);
    
    StartInterval=StartInterval(EmptyCheck);
    
    EndInterval=EndInterval(EmptyCheck);
    
    if sum(isnan(StartInterval))==0 && sum(isnan(EndInterval))==0
    tmp=cell(length( TrialIndex_tmp),1);
    try
    tmp(EmptyCheck)=cellfun(@(x,y) x(StartInterval((y)):min(EndInterval((y)),size(x,2))),EyeData(EmptyCheck),SelectionMat,'UniformOutput',0);
    catch
        keyboard
    end
     
    EyeInterval{i}=tmp;
       
    
    elseif sum(isnan(StartInterval))~=length(StartInterval) || sum(isnan(EndInterval))~=length(EndInterval)
        
        temp=cellfun(@(x,y) x(StartInterval(y):EndInterval(y)),EyeData(~isnan(StartInterval)&~isnan(EndInterval),:),SelectionMat(~isnan(StartInterval)&~isnan(EndInterval),:),'UniformOutput',0);
        EyeInterval_temp([~isnan(StartInterval)&~isnan(EndInterval)],:)=temp;
        EyeInterval{i}=EyeInterval_temp;
        
    else
        
        EyeInterval{i}=EyeInterval_temp;
       % disp('Fail to select interval,Empty!')
        
        
    end
    
    
    
end




end