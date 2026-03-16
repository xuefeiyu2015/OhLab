function TargetType=FractalTypeIdentifyBonus(TargetIndex)
%To identify whether an object is good or bad 
%First screen out the novel object, type III
%Temporally make it hardwired
TargetType=NaN*ones(size(TargetIndex));

NovelFractalStart=80;
TargetType(TargetIndex>=NovelFractalStart & TargetIndex<93)=3;%NovelObject
TargetType(TargetIndex>=13 & TargetIndex<=24)=1;%GoodObject
TargetType(TargetIndex>=25 & TargetIndex<=36)=2;%BadObject

TargetType(isnan(TargetType)&(mod(TargetIndex,2)==1))=1;%Good Object
TargetType(isnan(TargetType)&(mod(TargetIndex,2)==0))=2;%Bad Object
end