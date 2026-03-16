function SegmentData=EventSegTool(StartEvent,EndEvent,PlotStartOffset,PlotStopOffset_ori,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth)
%Function to segregate spike time according to the event marker, 
%the output of this function contains the mean firing rate beetween
%different events, the average PSTH for each pair of event bin, and the 
%PSTH for each pair of event bin for each trial,and the relative marker
%time bin;
%The output is:
%{
SegmentData.MeanFR=[AverageFiringRate_Event_Mean',AverageFiringRate_Event_Sem'];
SegmentData.MeanFRRaw

SegmentData.MeanPSTH=[PSTH_Event_Time',PSTH_Event'];
SegmentData.RelativeTimeMarker
SegmentData.AbsoluteTimeMarker

SegmentData.MeanFREach
SegmentData.MeanPSTHEach=[PSTH_Event_Each_Time',PSTH_Event_Each'];
SegmentData.RelativeTimeMarkerEach
SegmentData.AbosoluteTimeMarkerEach
%}




%by Xuefei Yu;
%04042020
%Intilization
ECode;
SegmentData=[];

PSTH_Event_Each=[];
PSTH_Event=[];
PSTH_Event_Time=[];
PSTH_Event_Each_Time=[];

EventTimeMarkerEach=[];
EventTimeMarker=[];
AbosoluteTimeMarkerEach=[];
AbosoluteTimeMarker=[];


AverageFiringRate_Event_Mean=[];
AverageFiringRate_Event_Sem=[];
AverageFiringRate_Event_Each=[];

RasterAfterAlign=[];
RasterAfterAlign_Each=[];
PlotStopOffset=[];


%Align Data with each events
for i=1:length(StartEvent)
    
   %Align with the current event
   
   AlignTime=FindOutTime(EventChannel,EventTimeChannel,StartEvent(i));
   
  
   
   
   %Initialize for the mean firing rate
   AverageFiringRateInInterval_Whole=[];
  %AverageFiringRateInInterval=[];
  SpikeFiringRate_temp=[];
  SpikeCount_temp=[];
  TimeSequence_temp=[];
  EventTime=[];
  AbosoluteTime=[];
 if ~sum(~isnan(AlignTime))
       disp('Current Event Marker Not Exist');
      AverageFiringRate_Event_Mean(i)=NaN;
      
      AverageFiringRate_Event_Sem(i)=NaN;
      AverageFiringRate_Event{i}=NaN;

PSTH_Event_Time{i}=NaN;
PSTH_Event{i}=NaN;
EventTimeMarker{i}=NaN;
AbosoluteTimeMarker{i}=NaN;

AverageFiringRate_Event_Each{i}{1}=NaN;
PSTH_Event_Each_Time{i}{1}=NaN;
PSTH_Event_Each{i}{1}=NaN;
EventTimeMarkerEach{i}{1}=NaN;
AbosoluteTimeMarkerEach{i}{1}=NaN;

RasterAfterAlign_Each{i}{1}=NaN;
RasterAfterAlign{i}=NaN;

PSTH_Event_SpikeCount{i}=NaN;
PSTH_Event_SpikeCount_Each{i}{1}=NaN;
       
 else
 
  for j=1:size( AlignTime,2)
       AlignTime_Curr=AlignTime(:,j);
       
       if isnan(EndEvent(i))
           if ~isnan(PlotStopOffset_ori(i))
              EndTime_Curr= AlignTime_Curr+PlotStopOffset_ori(i);
              PlotStopOffset(i)=PlotStopOffset_ori(i)/2;
           else
              EndTime_Curr= SpikeTime(:,end);
             PlotStopOffset(i)=0;
           end
    
       else
       EndTime=FindOutTime(EventChannel,EventTimeChannel,EndEvent(i));
       EndTime_Curr=EndTime(:,j);
       PlotStopOffset(i)=PlotStopOffset_ori(i);
    
       end
       
   
     %  [mindiff,Index]=min(abs(EndTime-AlignTime_Curr),[],2,'omitnan');
     %  EndTime_Curr=EndTime(Index);
     %{
     %For Robin081120FORAGEXPLO4.rst
     if(StartEvent(i)==SCENE_ON)
         %There shouldn't be an nan
         missingtrial=sum(isnan(AlignTime_Curr));
         if sum(isnan(AlignTime_Curr)) 
             if isnan(EndTime_Curr(isnan(AlignTime_Curr)))
                 normS=mode(AlignTime_Curr(~isnan(AlignTime_Curr)));
                 normE=mode(EndTime_Curr(~isnan(AlignTime_Curr)));
                 
                  AlignTime_Curr(isnan(AlignTime_Curr))=normS;
                  EndTime_Curr(isnan(AlignTime_Curr))=normE;
             else
                 AlignTime_Curr(isnan(AlignTime_Curr))=EndTime_Curr(isnan(AlignTime_Curr))-1000;%1000 is the imperical value for the scene onset from the fix onset
     
             end
             
             fprintf('%d scene start code is missing, the current trial scene onset maybe inaccurate\n',missingtrial);
         end
         
     end
    %}
       
       %Remove the nans
       notanflag=~isnan(AlignTime_Curr);
       AlignTime_Curr=AlignTime_Curr(notanflag);
       EndTime_Curr=EndTime_Curr(notanflag);
       
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
  
       
       %Select interval to compute the average firng rate 
       AverageFRTime=[AlignTime_Curr-AlignTime_Curr,EndTime_Curr-AlignTime_Curr];%in ms;%Time interval to compute average firing rate within the time bin
       if isnan(PlotStartOffset(i))
           StartOffset=SpikeTime(:,end-1)-AlignTime_Curr;
       else
           StartOffset=-PlotStartOffset(i)*ones(size(AlignTime_Curr,1),1);
           
       end
       
       %Time window to plot the PSTH      
       SelectIntervalTime=[StartOffset,EndTime_Curr-AlignTime_Curr+PlotStopOffset(i)*ones(size(AlignTime_Curr,1),1)];
       
       %Align spikes according to the marker
       SpikeTimeAfterAlign_All=SpikeTime(notanflag,:)-AlignTime_Curr;
      % re_ISI_All = [1./diff(SpikeTimeAfterAlign_All,1,2),0*ones(size(SpikeTimeAfterAlign_All,1),1)]*1000; %in Hzs
       [re_ISI{j},re_ISI_seq{j}] = ISIInInterval(SelectIntervalTime,SpikeTimeAfterAlign_All,BinWidth);
       
       
         
       %Using a gaussian smooth window to calculate firing rate
       
        [SpikeFiringRate_temp{j},TimeSequence_temp{j},SpikeCount_temp{j}]=CalculateSpikeFiringRate(SpikeTimeAfterAlign_All,StepSize,BinWidth,SelectIntervalTime);
       
        %Select spike time interval along with the PSTH
       SpikeTimeInSelectInterval=SelectSpikeTimeInterval(SpikeTimeAfterAlign_All,SelectIntervalTime(:,1),SelectIntervalTime(:,2));%Return spike count to calculate average firing rate
      
      %Measure the average firing rate
       SpikeTimeInInterval=SelectSpikeTimeInterval(SpikeTimeAfterAlign_All,AverageFRTime(:,1),AverageFRTime(:,2));%Return spike count to calculate average firing rate
       SpikeCountInInterval=sum(~isnan(SpikeTimeInInterval),2);
       
       %Check whether out of interval
       
       SpikeTimeStart=SpikeTimeAfterAlign_All(:,end-1);
       SpikeTimeEnd=SpikeTimeAfterAlign_All(:,end);
       SpikeRateBin_left=max([SpikeTimeStart,AverageFRTime(:,1)],[],2);
       SpikeRateBin_right=min([SpikeTimeEnd,AverageFRTime(:,2)],[],2);
       
       
       SpikeRateInInterval=SpikeCountInInterval./((SpikeRateBin_right-SpikeRateBin_left)/1000);
       
       AverageFiringRateInInterval_Whole=[AverageFiringRateInInterval_Whole;SpikeRateInInterval];
       
       %Store the average firing rate and PSTH
      
       AverageFiringRate_Event_Each{i}{j}=SpikeRateInInterval;%in Hz
      

       PSTH_Event_Each{i}{j}=SpikeFiringRate_temp{j};
       PSTH_Event_Each_Time{i}{j}=TimeSequence_temp{j};
       %Store the PSTH spike count
       PSTH_Event_SpikeCount_Each{i}{j}=SpikeCount_temp{j};
       
       %Store the relative time stamp
        EventTime=[EventTime;AverageFRTime];
       
       EventTimeMarkerEach{i}{j}=AverageFRTime;
       
       %Store the abosolute time stamp
       
       AbosoluteTime=[AbosoluteTime;AverageFRTime+AlignTime_Curr];
       AbosoluteTimeMarkerEach{i}{j}=AbosoluteTime;    
       %Store the raster
       RasterAfterAlign_Each{i}{j}= SpikeTimeInSelectInterval;
       
       SpikeTimeInSelectInterval_temp{j}= SpikeTimeInSelectInterval;
    
  end
 
   %Store the PSTH
   temp_length_max=max(cellfun(@(x) size(x,2),SpikeFiringRate_temp));
   
   PSTH_Event{i}=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),temp_length_max-size(x,2))],SpikeFiringRate_temp,'uniform',0)');
   PSTH_Event_Time{i}=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),temp_length_max-size(x,2))],TimeSequence_temp,'uniform',0)');
   EventTimeMarker{i}=EventTime;
   AbosoluteTimeMarker{i}=AbosoluteTime;
   %Store the spike count PSTH
   
   PSTH_Event_SpikeCount{i}=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),temp_length_max-size(x,2))],SpikeCount_temp,'uniform',0)');
  
   
    %Store the average firing rate   
    AverageFiringRate_Event{i}=AverageFiringRateInInterval_Whole;

    AverageFiringRate_Event_Mean(i)=nanmean(AverageFiringRateInInterval_Whole);
    AverageFiringRate_Event_Sem(i)=nanstd(AverageFiringRateInInterval_Whole)/sqrt(length(AverageFiringRateInInterval_Whole(~isnan(AverageFiringRateInInterval_Whole))));
    
    %Store the raster
    temp_length_max=max(cellfun(@(x) size(x,2), SpikeTimeInSelectInterval_temp));
   
     RasterAfterAlign{i}=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),temp_length_max-size(x,2))], SpikeTimeInSelectInterval_temp,'uniform',0)');
  end %End of the marker existence check
end

%Organize the output

%SegmentData.
SegmentData.MeanFR=[AverageFiringRate_Event_Mean',AverageFiringRate_Event_Sem'];
SegmentData.MeanFRRaw=AverageFiringRate_Event;

SegmentData.reISI = re_ISI;
SegmentData.reISIseq = re_ISI_seq;





SegmentData.MeanPSTH=[PSTH_Event_Time',PSTH_Event'];
SegmentData.RelativeTimeMarker=EventTimeMarker;
SegmentData.AbsoluteTimeMarker=AbosoluteTimeMarker;

SegmentData.MeanFREach=AverageFiringRate_Event_Each;
SegmentData.MeanPSTHEach=[PSTH_Event_Each_Time',PSTH_Event_Each'];
SegmentData.RelativeTimeMarkerEach=EventTimeMarkerEach;
SegmentData.AbsoluteTimeMarkerEach=AbosoluteTimeMarkerEach;

SegmentData. RasterAfterAlignEach= RasterAfterAlign_Each;
SegmentData. RasterAfterAlign= RasterAfterAlign;

SegmentData.PSTH_Event_SpikeCount=PSTH_Event_SpikeCount';
SegmentData.PSTH_Event_SpikeCountEach=PSTH_Event_SpikeCount_Each;


end