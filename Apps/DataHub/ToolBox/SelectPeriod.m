function Select=SelectPeriod(Original,Interval,ToleranceMode);
%Function to select the data within the interval
%Xuefei Yu 04212020
for i=1:length(Original)
%for i=252:252
    data_curr=Original{i};
    if isempty(data_curr)
        %Select{i}=data_curr;
        Select{i}=[];
    else
    Interval_curr=Interval(i,:);
    Selection=find((data_curr(1,:)>=Interval_curr(1)) &(data_curr(1,:)<=Interval_curr(2)));
    if isempty(Selection)
        
        Selection_temp=find(Interval_curr(1)>=data_curr(1,:),1,'first');
        if isempty(Selection_temp) || Selection_temp==1
            if ToleranceMode
            Selection=1;
            else
            Selection=[];
            end
            
        else
            if ToleranceMode
        
             Selection=[Selection_temp,Selection_temp+1];
     
            else
                Selection=[];
            end
                
        end
    end
    tmp=data_curr(:,[Selection]);
    if ToleranceMode==2 %Strick interval
        tmp(1,1)=max(tmp(1,1),Interval_curr(1));
        tmp(end,end)=min(tmp(end,end),Interval_curr(2));
    end
     Select{i}=tmp;
    end
    
end


end