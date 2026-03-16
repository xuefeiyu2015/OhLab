function currContainIndexCorrect=BonusTrialConditionScreen(TargetIndex,TargetCondition,ContainCoupleIndex);
%This function is to screen out the extra trial conditions according to the
%target 
%Xuefei Yu 02192020
Index=TargetIndex;
currContainIndexCorrect=[];
for i=1:length(Index)
   currIndex=Index(i);
   currContainIndex=ContainCoupleIndex(ismember(ContainCoupleIndex(:,1),currIndex),:);
   currTargetIndex=TargetCondition(TargetCondition==currIndex,:);
    if(size(currContainIndex,1)~=size( currTargetIndex,1))
       currContainIndexCorrect=[currContainIndexCorrect;currContainIndex(size(currContainIndex,1)-size( currTargetIndex,1)+1:end,:)];
   else
       
       currContainIndexCorrect=[currContainIndexCorrect;currContainIndex];
       
       
   end
    
end

end