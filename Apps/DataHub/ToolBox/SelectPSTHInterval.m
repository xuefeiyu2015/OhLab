function [SelctedPSTH,SelctedPSTH_Time]=SelectPSTHInterval(PSTH,PSTH_time_ori,Interval_ori)
%Function to select PSTH intervals;

SelctedPSTH=NaN*ones(size(PSTH));
SelctedPSTH_Time=NaN*ones(size(PSTH));
MaxLength=0;

PSTH_time=PSTH_time_ori-Interval_ori(1);
Interval=Interval_ori-Interval_ori(1);

Check=PSTH_time(:,1)>Interval(1);
if sum(Check)>0
    MaxFill=max(PSTH_time(Check,1));
    diff_TimeBin=nanmean((PSTH_time_ori(:,2)-PSTH_time_ori(:,1)));
    numOfFill=ceil(MaxFill/diff_TimeBin)+1;
    
    for i=1:size(PSTH_time,1)
        if isnan(PSTH_time(i,1))
            FillTimeBin(i,:)=NaN*ones(numOfFill,1);
        else
            
            FillTimeBin(i,:)=PSTH_time(i,1)-numOfFill*diff_TimeBin:diff_TimeBin:PSTH_time(i,1)-diff_TimeBin;
        end
        
        
    end
    
    PSTH_time=[FillTimeBin,PSTH_time];
    PSTH=[NaN*FillTimeBin,PSTH];
    
    
   % fprintf('Auto fill in by %d time bins with %d time interval \n',numOfFill,diff_TimeBin);
end


for i=1:size(PSTH,1)
   Select=PSTH_time(i,:)<=Interval(2) & PSTH_time(i,:)>=Interval(1);
   if ~sum(Select==1)
       %If nothing in the interval,use the mean of the others
      Select=nanmean(PSTH_time,1)<=Interval(2) & nanmean(PSTH_time,1)>=Interval(1); 
       
   end
   
   SelctedPSTH(i,1:sum(Select))=PSTH(i,Select);
    SelctedPSTH_Time(i,1:sum(Select))=PSTH_time(i,Select);
    if MaxLength<sum(Select)
        MaxLength=sum(Select);
    end
    
end

SelctedPSTH=SelctedPSTH(:,1:MaxLength);
SelctedPSTH_Time=SelctedPSTH_Time(:,1:MaxLength)+Interval_ori(1);

end