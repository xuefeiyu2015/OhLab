function AbnormalTrials=AbnormalTrialCheck(Data,Marker);
%This is the abnormal trial check, if it returns sth, there is abnormal
%trials mainly due to the photodio error; such trials should be removed 
ECode;
EventChannel=Data.EventChannel;
EventTimeChannel=Data.EventTimeChannel;
if length(Marker)==1
    %Single marker check
CheckTime=FindOutTime(EventChannel,EventTimeChannel,Marker);

AbnormalTrials=find(isnan(CheckTime)==1);
if ~isempty(AbnormalTrials)
    fprintf('Warning!!Abnormal trial detected!/n Trial %d will be removed!\n',AbnormalTrials);
end
elseif length(Marker)==2
    %Marker Pair Check
    CheckTimePair1=~isnan(FindOutTime(EventChannel,EventTimeChannel,Marker(1)));
    CheckTimePair2=~isnan(FindOutTime(EventChannel,EventTimeChannel,Marker(2)));

    if size(CheckTimePair2,2) > size(CheckTimePair1,2)
        CheckTimePair2 = CheckTimePair2(:,1);
    end
    
    
     AbnormalTrials=find(sum(~(CheckTimePair1 == CheckTimePair2),2)==1);
     
   if ~isempty(AbnormalTrials)
    fprintf('Warning!!Abnormal trial detected!/n Trial %d will be removed!\n',AbnormalTrials);
   end
else
    %Equality check
    AbnormalTrials=[];
    for i = 1:size(Marker,2)
        Pair = [Marker(1,i),Marker(2,i)];
        CheckTimePair1=~isnan(FindOutTime(EventChannel,EventTimeChannel,Pair(1)));
        CheckTimePair2=FindOutTime(EventChannel,EventTimeChannel,Pair(2));
        if size(CheckTimePair2,2) > size(CheckTimePair1,2)
            
            CheckTimePair2 = ~isnan(CheckTimePair2(:,i));
            
        else
            CheckTimePair2 = ~isnan(CheckTimePair2);
            
            
            
        end
        
         AbnormalTrials=[AbnormalTrials,find(sum(~(CheckTimePair1 == CheckTimePair2),2)==1)];
     
  
        
    end
     if ~isempty(AbnormalTrials)
        fprintf('Warning!!Abnormal trial detected!/n Trial %d will be removed!\n',AbnormalTrials);
    end
        
        
    
end 
    
end
