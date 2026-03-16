 function Outcome=TransformToOutcome(TargetIndex);
 %Function to transform outcome intp rejection and acceptions ;
 %1:Accept;0:Reject
 %Xuefei Yu 05022020
  Outcome=0*ones(size(TargetIndex));
   for ii=1:size(TargetIndex,1)
      EndPoint=find(isnan(TargetIndex(ii,:))==1,1,'first')-1;
      Outcome(ii,EndPoint)=1;%Acception
      
      Outcome(ii,isnan(TargetIndex(ii,:)))=NaN;
      
            
            
            
   end
end