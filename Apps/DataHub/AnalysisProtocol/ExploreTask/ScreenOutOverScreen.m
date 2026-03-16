function FixOffTime_screened=ScreenOutOverScreen(FixOffTime,OutScreenTime);
%Function to screen out the fix off time for overscreen trials 

OutScreenFlag=~isnan(OutScreenTime);
FixOffTime_OverScreen=FixOffTime(OutScreenFlag,:);
for i=1:size(FixOffTime_OverScreen,1)
   tmp=FixOffTime_OverScreen(i,:);
   idx=find(~isnan(tmp)==1);
   FixOffTime_OverScreen(i,idx(end))=NaN;
end


FixOffTime_screened=FixOffTime;
FixOffTime_screened(OutScreenFlag,:)=FixOffTime_OverScreen;
%Remove the whole nan col
Col_NaN=find(sum(isnan(FixOffTime_screened),1)==size(FixOffTime_screened,1));
FixOffTime_screened(:,Col_NaN)=[];

end