function StartTime = FindOutTime(Event,EventTime,Marker)
% This function is to find out the event time according to the marker name
%Event is the event channel; eventtime is the event time channel; marker is
%the marker name of interest
%It can be used to findout anything in the EventTime channel according to
%the marker in the Event;
%Xuefei Yu 01032020
MaxLength=1;
StartTime=NaN*ones(size(Event));
for i=1:size(Event,1)
    
    Bin=find(Event(i,:)==Marker);
    
    StartTime(i,1:length(Bin))=EventTime(i,Bin);
    if(length(Bin)>MaxLength)
        MaxLength=length(Bin);
    end
    
    
    
end

StartTime=StartTime(:,1:MaxLength);

end

