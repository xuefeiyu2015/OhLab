function [FiringRate,TimeSequence,SpikeCount]= CalculateSpikeFiringRate(SpikeTime_Ori,StepSize,BinWidth_Smooth,TimeInterval)
%Transform Spike Time into Spike Firing Rate; YXF 10152019
%   With a sliding window
SpikeTime=SpikeTime_Ori(:,1:end-2);
SpikeTimeEnd=SpikeTime_Ori(:,end);%Spike End time

SpikeTimeStart=SpikeTime_Ori(:,end-1);%Spike Start time

    



MaxLength=0;
NumForSmoothOffset=5;
BinWidth=StepSize;%Using gaussian smooth to control the overlap

FiringRate=NaN*ones(size(SpikeTime,1),ceil(nanmax(TimeInterval(:,2)-TimeInterval(:,1))/StepSize)+10);

TimeSequence=FiringRate;
SpikeCount=FiringRate;





for i=1:size(SpikeTime,1)
    
    %In each trial
    currentSpikeTime=SpikeTime(i,:);
  %Initialization
   SelectInterval=[];
   SelectTimeInterval=[];
   SelectIntervalBackward=[];
   SelectIntervalForward=[];
   TimeBackward=[];
   TimeForward=[];
   TimeIntervalShift=0;
   
   if TimeInterval(i,1)>0
       TimeIntervalShift=TimeInterval(i,1);     
       TimeInterval(i,2)=TimeInterval(i,2)-TimeInterval(i,1);
       TimeInterval(i,1)=0;
      
       
   end
   
    currentSpikeTime=currentSpikeTime-TimeIntervalShift;
    SpikeTimeStart(i)=SpikeTimeStart(i)-TimeIntervalShift;
    SpikeTimeEnd(i)=SpikeTimeEnd(i)-TimeIntervalShift;
   
      SmoothOffset=BinWidth* NumForSmoothOffset;%in ms; for smooth
      
      TimeInterval(i,1)=TimeInterval(i,1)-SmoothOffset;
      TimeInterval(i,2)=TimeInterval(i,2)+SmoothOffset;
      

    if TimeInterval(i,2)>0
       StepCountForward=max(1,ceil((TimeInterval(i,2)-BinWidth)/StepSize));
    else
        StepCountForward=ceil((TimeInterval(i,2)-BinWidth)/StepSize);
        
    end
    StepCountBackward=ceil((0-TimeInterval(i,1)-BinWidth)/StepSize);
    if MaxLength<((StepCountForward+StepCountBackward)+2-NumForSmoothOffset*2)
        MaxLength=(StepCountForward+StepCountBackward)+2-NumForSmoothOffset*2;
        
    
    end
   
    %First forward 
    if StepCountForward>0
        
      
        
        for j=0:StepCountForward
            if j*StepSize+BinWidth>TimeInterval(i,2)
                if j*StepSize+BinWidth>SpikeTimeEnd(i)
                    SelectIntervalForward(j+1)=NaN;
                else
                
                   SelectIntervalForward(j+1)=sum((currentSpikeTime<TimeInterval(i,2))&(currentSpikeTime>(j*StepSize)))/(TimeInterval(i,2)-(j*StepSize));
               end
            else
                if j*StepSize+BinWidth>SpikeTimeEnd(i)
                    SelectIntervalForward(j+1)=NaN;
                else
                    SelectIntervalForward(j+1)=sum((currentSpikeTime<(j*StepSize+BinWidth))&(currentSpikeTime>(j*StepSize)))/BinWidth;
                end
               
            end
             TimeForward(j+1)=j*StepSize+StepSize/2;
     
        end
        SelectInterval=SelectIntervalForward;
      
        SelectTimeInterval=TimeForward;
   
   %Then backward
  
    if StepCountBackward>0
        for j=linspace(0,-StepCountBackward, StepCountBackward+1)
            if j*StepSize-BinWidth<TimeInterval(i,1)
                if j*StepSize-BinWidth<SpikeTimeStart(i)
                    SelectIntervalBackward(-j+1)=NaN;
                else
                   SelectIntervalBackward(-j+1)=sum((currentSpikeTime<(j*StepSize))&(currentSpikeTime>TimeInterval(i,1)))/(-TimeInterval(i,1)-(j*StepSize));
                end
            else
                if j*StepSize-BinWidth<SpikeTimeStart(i)
                     SelectIntervalBackward(-j+1)=NaN;                    
                else               
                    SelectIntervalBackward(-j+1)=sum((currentSpikeTime<(j*StepSize))&(currentSpikeTime>(j*StepSize-BinWidth)))/BinWidth;
                end
               
            end
                
              TimeBackward(-j+1)=j*StepSize-StepSize/2;
 
        end
       
     %  SelectInterval_temp=SelectInterval;
         
         SelectInterval=[fliplr(SelectIntervalBackward),SelectInterval];
     %     SelectTimeInterval_temp=SelectTimeInterval;
         SelectTimeInterval=[fliplr(TimeBackward),SelectTimeInterval];
       
    end
   
    end
  %  SelectInterval_Smooth=[ smoothdata(fliplr(SelectIntervalBackward),'gaussian',BinWidth_Smooth/StepSize,'omitnan') smoothdata(SelectInterval,'gaussian',BinWidth_Smooth/StepSize,'omitnan')];
     SpikeCount(i,1:length(SelectInterval)-NumForSmoothOffset*2)= ceil(SelectInterval(NumForSmoothOffset+1:end-NumForSmoothOffset)*BinWidth);
     
    if ~isempty(SelectInterval)
        

       tmp=GaussSmooth(SelectTimeInterval,SelectInterval,BinWidth_Smooth)*1000;
       
       %Remove the smooth offset
           
       FiringRate(i,1:length(SelectInterval)-NumForSmoothOffset*2)=tmp(NumForSmoothOffset+1:end-NumForSmoothOffset);
       
       
    end
   
%    FiringRate(i,1:length(SelectInterval))= smoothdata(SelectInterval,'gaussian',BinWidth_Smooth/StepSize,'omitnan')*1000;

    TimeSequence(i,1:length(SelectInterval)-NumForSmoothOffset*2)=SelectTimeInterval(NumForSmoothOffset+1:end-NumForSmoothOffset)+TimeIntervalShift;

 
%   end %End if the isnan 
    
end
FiringRate=FiringRate(:,1:MaxLength);
TimeSequence= TimeSequence(:,1:MaxLength);
SpikeCount=SpikeCount(:,1:MaxLength);

end

