function alignedEye = AlignEyeDataWithEvents(StartEvent,EndEvent,PlotStartOffset,PlotStopOffset_ori,EventChannel,EventTimeChannel,EyeData);

%Reorganize eye data
Eye_X = ReorganizeEye(EyeData.X);
Eye_Y = ReorganizeEye(EyeData.Y);
Eye_Bin = EyeData.Bin;

EyeAligned_X ={};
EyeAligned_Y ={};
EyeTimeAligned={};

%Align Data with each events
for i=1:length(StartEvent)
    

 AlignTime=FindOutTime(EventChannel,EventTimeChannel,StartEvent(i));
 if ~sum(~isnan(AlignTime))
       disp('Current Event Marker Not Exist');
      
 else
     for j=1:size( AlignTime,2)
         AlignTime_Curr=AlignTime(:,j);
         
        if isnan(EndEvent(i))
           if ~isnan(PlotStopOffset_ori(i))
              EndTime_Curr= AlignTime_Curr+PlotStopOffset_ori(i);
              PlotStopOffset(i)=PlotStopOffset_ori(i)/2;
           else
              EndTime_Curr= EyeData.MaxTime;
             PlotStopOffset(i)=0;
           end
    
       else
        EndTime=FindOutTime(EventChannel,EventTimeChannel,EndEvent(i));
        EndTime_Curr=EndTime(:,j);
        PlotStopOffset(i)=PlotStopOffset_ori(i);  
        end
       %{
        %Remove the nans
       notanflag=~isnan(AlignTime_Curr);
       AlignTime_Curr=AlignTime_Curr(notanflag);
       EndTime_Curr=EndTime_Curr(notanflag);
       %}
       %Check whether end time is later than the start time
       
       Check=EndTime_Curr-AlignTime_Curr;
       NegativeTime=Check<0;
       
       if sum(NegativeTime)>0
           MessageStr=sprintf('Correct negative duration with the mid value for event %d, term %d',StartEvent(i),j);
           disp(MessageStr);
           CorrectNum=sum(NegativeTime);
           MessageStr2=sprintf('Number of correction: %d',CorrectNum);
           
       
       AllNormalTime=Check(~NegativeTime);
       AllNormalTime_Sort=sort(AllNormalTime);
       
       Mid_AllNormalTime_Sort=AllNormalTime_Sort(floor(length(AllNormalTime_Sort)/2));
       
       
       EndTime_Curr(NegativeTime)=AlignTime_Curr(NegativeTime)+Mid_AllNormalTime_Sort;
       end
       %Align eye data
       EyeTime = 1:Eye_Bin:size(Eye_X,2);
       EyeAligned_X_tmp={};
       EyeAligned_Y_tmp={};
       TimeAfterAlign=[];
       
       for m = 1:size(Eye_X,1)
           
           TimeAfterAlign(m,:)=EyeTime-AlignTime_Curr(m);
           
           Interval =[-PlotStartOffset(i), EndTime_Curr(m)-AlignTime_Curr(m)+PlotStopOffset(i)];
           
           if TimeAfterAlign(m,1) >Interval(1)
               select_start = max(EyeTime(1)-AlignTime_Curr);
               sel = TimeAfterAlign(m,:)>=select_start & TimeAfterAlign(m,:)<=Interval(2);

           else
               sel = TimeAfterAlign(m,:)>=Interval(1) & TimeAfterAlign(m,:)<=Interval(2);
               
           end
           
           
           
           EyeAligned_X_tmp{m} = Eye_X(m,sel);
           EyeAligned_Y_tmp{m} = Eye_Y(m,sel);
          
           
           TimeAfterAlign_tmp{m} = TimeAfterAlign(m,sel);
           
           
       end
       if size( AlignTime,2)==1
           EyeAligned_X{i}=ReorganizeEye(EyeAligned_X_tmp);
           EyeAligned_Y{i}=ReorganizeEye(EyeAligned_Y_tmp);
           EyeTimeAligned{i} =ReorganizeEye(TimeAfterAlign_tmp);
       else
        EyeAligned_X{i,j}=ReorganizeEye(EyeAligned_X_tmp);
        EyeAligned_Y{i,j}=ReorganizeEye(EyeAligned_Y_tmp);
        EyeTimeAligned{i,j} =ReorganizeEye(TimeAfterAlign_tmp);
       end
       
       
       
       
       
         
         
         
     end
     
 end
 
 
 %AlignedData=AlignEyeData(Data,EyeTime,MarkerTime,Interval);
   
end



alignedEye.X=EyeAligned_X;
alignedEye.Y=EyeAligned_Y;
alignedEye.Time = EyeTimeAligned;



end
