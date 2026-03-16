function OutputTime=FindOutLatency(Time,Cri,Threshold)
%Find out the time according to the criterium,
%Xuefei Yu
%04052020

BinSatisfied=find(Cri==1)';

if Threshold==1
    OutputTime=Time(BinSatisfied(1));
else
    if length(BinSatisfied)<Threshold
        OutputTime=[];
    else
   Mat=repmat(BinSatisfied,1,Threshold);
   for i=1:Threshold
       MatTrans(:,i)=[Mat(i:length(Mat)-Threshold+i,i)];
   end
   
   OutputTime=Time(BinSatisfied(find(sum(diff(MatTrans,1,2)==1,2)==Threshold-1,1,'first')));
    end
   
end





end