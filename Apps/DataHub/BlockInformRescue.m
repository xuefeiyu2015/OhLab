  %Rescue block information if there is some error in data saving 
  %Xuefei Yu, 02052020
function BlockType=BlockInformRescue(BlockType_Unrescued,LostTrials,FixationPointIndex,SceneAndObject);
ECode;

FixationPointIndexLost=FixationPointIndex(LostTrials,2);

BlockGuess=~(FixationPointIndexLost>2);


SceneAndObjectLost=SceneAndObject(ismember(SceneAndObject(:,1),LostTrials),:);
UniqueScene=unique(SceneAndObjectLost(:,3));
if length(UniqueScene)~=6
    SceneType=999*ones(length(LostTrials),1);
else
    SceneType=0*ones(length(LostTrials),1);
    SceneType((FixationPointIndexLost==FIXFRB),:)=1;
    

end



BlockRescue=[LostTrials,BlockGuess,SceneType];
BlockType=[BlockRescue;BlockType_Unrescued;];
end