function Batch_FourBlockCmp_Exploration(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
%Batch analysis for comparison of the scene selectivity index for all the
%four blocks; by Xuefei Yu 03172021
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesPath=Data.ResultFilePath(StartFile: EndFile);
FilesName=Data.FileName;
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
TotalNeuronNum=sum(cellfun(@numel,FilesPath));
ChannelID=Data.ChannelNumber;

PlotExample=1;
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
         
StableSceneColor=[ 1.00 0.54 0.00];
FlexibleSceneColor=[ 0.25 0.80 0.54];
DifferenceColor=[0.4471    0.3020    0.9176];

ChannelFull=[];


%Load file
FileIndex=1;
for i=1:length(FilesPath)
  clear OutputData;
  ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];
  
  for j=1:length(FilesPath{i})
      currFileName=FilesPath{i}{j};
      load(currFileName);
     OutputDataTemp=OutputData; 
     
     Task(FileIndex)=OutputData.FourBlockSceneCompare(1).TaskCode;
     
     TaskName{FileIndex}=OutputData.FourBlockSceneCompare(1).Task;  
     FilesNameEach{FileIndex}=FilesName{i};
     
    
     %Separate Keyvalues in each container;
    % Stamp_tmp=OutputData.EventSegRT(1).DataStamp;
    % KeyWhole=keys(Stamp_tmp);
     
     FourBlockSceneCompare{FileIndex}=OutputData.FourBlockSceneCompare(1).DataStamp;
     
      FileIndex=FileIndex+1;
      
  end
  
    
     
     
end
%Assume that keys are all the same for each file, use the first file to
%abstract keys for all files;
KeysLib=keys(FourBlockSceneCompare{1});
AllValues=cellfun(@(x) values(x),FourBlockSceneCompare,'uniform',0);
%Reorganize the value according to each key
for i=1:length(KeysLib)
  %  currKey=KeysLib{i};
   KeyValues{i}=cellfun(@(x) x{i},AllValues,'uniform',0);
 
end

SSISceneKey={'Time_PSTH_SceneOn','Time_PSTH_SceneOff','SSI_Group','SSI_BlockIndex','SceneTypeIndex','SSI_GroupOff',' p_SSI_GroupOff',' p_SSI_GroupOff',... 
    };

%Time 
SSISceneOnTime=KeyValues{ismember(KeysLib,SSISceneKey{1})};
 
 SSISceneOn_TimeMat=fillinnancol(SSISceneOnTime);
 SSISceneOn_TimeMean=nanmean(SSISceneOn_TimeMat,1);
 
 SSISceneOffTime=KeyValues{ismember(KeysLib,SSISceneKey{2})};
 
 SSISceneOff_TimeMat=fillinnancol(SSISceneOffTime);
 SSISceneOff_TimeMean=nanmean(SSISceneOff_TimeMat,1);
 
%SSI
SceneGroup={[1000:1007],[1010:1017],[1054:1057,1060:1063]};
SceneGroupStr={'FlexibleScenes','StableScenes','DuelScenes'};

SSIOnData=KeyValues{ismember(KeysLib,SSISceneKey{3})};
SSIOnBlockData=KeyValues{ismember(KeysLib,SSISceneKey{4})};
SSISceneTypeData=KeyValues{ismember(KeysLib,SSISceneKey{5})};
%SSI_SceneOnTime=cellfun(@(x) cellfun(@(y) y(:,[2:3]),x,'uniform',0),KeyValues('Time_PSTH_SceneOn'),'uniform',0);%Normalize
% PSTHSceneData_Time=cellfun(@(x) cellfun(@(y) y(:,1),x,'uniform',0),KeyValues(Sel_tmp),'uniform',0);
 



%First compare the stable scene with the flexible scene with the dual scene
%for both fix location and random location conditions
SelectionScene=cellfun(@(x) sum(x)==3,SSISceneTypeData);
sel_tmp=cellfun(@(x) SelectLocCompare(x),SSIOnBlockData);
SelectionBlockType=arrayfun(@(x) x.flag,sel_tmp);
SelectionBlockIndex=arrayfun(@(x) x.index,sel_tmp,'uniform',0);


FixLocRandLocThreeCompare=SelectionScene & SelectionBlockType;

SSIOnD_LocTri=SSIOnData(FixLocRandLocThreeCompare);
BlockTri=SSIOnBlockData(FixLocRandLocThreeCompare);
SelBlockIndex=SelectionBlockIndex(FixLocRandLocThreeCompare);

SSIOnSelection_LocTri=cellfun(@(x,y) selectdata(x,y),SSIOnD_LocTri,SelBlockIndex,'uniform',0);
SSIOn_FleTri=cellfun(@(x) x{1}{1},SSIOnSelection_LocTri,'uniform',0);
SSIOn_FleTri_Mat=fillinnancol(SSIOn_FleTri);


SSIOn_StaTri=cellfun(@(x) x{2}{1},SSIOnSelection_LocTri,'uniform',0);
SSIOn_StaTri_Mat=fillinnancol(SSIOn_StaTri);

SSIOn_DuFleTri=cellfun(@(x) x{3}{1},SSIOnSelection_LocTri,'uniform',0);
SSIOn_DuFleTri_Mat=fillinnancol(SSIOn_DuFleTri);

SSIOn_DuStaTri=cellfun(@(x) x{3}{2},SSIOnSelection_LocTri,'uniform',0);
SSIOn_DuStaTri_Mat=fillinnancol(SSIOn_DuStaTri);


SSIOn_FleTri_Mat_Mean=nanmean(SSIOn_FleTri_Mat,1);
SSIOn_StaTri_Mat_Mean=nanmean(SSIOn_StaTri_Mat,1);

SSIOn_FleTri_Mat_Sem=nanstd(SSIOn_FleTri_Mat,[],1)./sqrt(sum(~isnan(SSIOn_FleTri_Mat),1));
SSIOn_StaTri_Mat_Sem=nanstd(SSIOn_StaTri_Mat,[],1)./sqrt(sum(~isnan(SSIOn_StaTri_Mat),1));


SSIOn_DuFleTri_Mat_Mean=nanmean(SSIOn_DuFleTri_Mat,1);
SSIOn_DuStaTri_Mat_Mean=nanmean(SSIOn_DuStaTri_Mat,1);

SSIOn_DuFleTri_Mat_Sem=nanstd(SSIOn_DuFleTri_Mat,[],1)./sqrt(sum(~isnan(SSIOn_DuFleTri_Mat),1));
SSIOn_DuStaTri_Mat_Sem=nanstd(SSIOn_DuStaTri_Mat,[],1)./sqrt(sum(~isnan(SSIOn_DuStaTri_Mat),1));


%Then compare the four blocks for duel scenes

SelectionSceneDuel=cellfun(@(x) x(3)==1,SSISceneTypeData);
sel_tmp=cellfun(@(x) SelectSeqCompare(x),SSIOnBlockData);
SelectionBlockTypeII=arrayfun(@(x) x.flag,sel_tmp);
SelectionBlockIndexII=arrayfun(@(x) x.index,sel_tmp,'uniform',0);

Duel4Compare=SelectionSceneDuel & SelectionBlockTypeII;

SSIOn_Duel4C=SSIOnData(Duel4Compare);
BlockDuel4=SSIOnBlockData(Duel4Compare);
SelBlockIndexDuel4=SelectionBlockIndexII(Duel4Compare);


SSIOnSelection_Seq=cellfun(@(x,y) selectdata(x,y),SSIOn_Duel4C,SelBlockIndexDuel4,'uniform',0);

SSIOn_Duel4=cellfun(@(x) x{3}(3),SSIOnSelection_Seq,'uniform',0);




FigureStartNum=300;
FigureIndex=1;
if ShowFigureFlag
    
    figtitlestr{FigureIndex}='SSI_Loc_StaFleDuel';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,800],'Name',figtitlestr{FigureIndex});
    subplot(2,2,1);
    errorbar(SSISceneOn_TimeMean,SSIOn_FleTri_Mat_Mean,SSIOn_FleTri_Mat_Sem,'-b');
    hold on
    errorbar(SSISceneOn_TimeMean,SSIOn_StaTri_Mat_Mean,SSIOn_StaTri_Mat_Sem,'-r');
    subplot(2,2,2);
    errorbar(SSISceneOn_TimeMean,SSIOn_DuFleTri_Mat_Mean,SSIOn_DuFleTri_Mat_Sem,'-b');
    hold on
    errorbar(SSISceneOn_TimeMean,SSIOn_DuStaTri_Mat_Mean,SSIOn_DuStaTri_Mat_Sem,'-r');
    subplot(2,2,3);
    imshow([SSIOn_StaTri_Mat,SSIOn_FleTri_Mat]);
    colormap jet;
    subplot(2,2,4);
    imshow([SSIOn_DuStaTri_Mat,SSIOn_DuFleTri_Mat]);
    colormap jet;
    
 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;   


end


keyboard

end
function mat=fillinnancol(mat_ori);
MaxCol=max(cellfun(@(x) size(x,2),mat_ori));
mat=cell2mat(cellfun(@(x) [x, NaN*ones(size(x,1),MaxCol-size(x,2))]',mat_ori,'uniform',0))';

end
function sel=SelectLocCompare(x);
%Select out conditions for fix loc and rand loc comparison
filter={[1],[11],[1,11]};
cri=zeros(1,3);
index=[];
for i=1:3
    x_tmp=cell2mat(x(i,:));
    filter_curr=filter{i};
    test=ismember(x_tmp,filter_curr);
    if sum(test)==length(filter_curr)
        cri(i)=1;
        index{i}=test;
    end
    
end

if sum(cri)==3
    sel.flag=1;
    sel.index=index;
else
    sel.flag=0;
    sel.index=[];
end


end
function sel=SelectSeqCompare(x);
%Select out conditions for fix loc and rand loc comparison
filter={[],[],[1,11,10,0]};
cri=zeros(1,3);
index=[];
for i=1:3
    x_tmp=cell2mat(x(i,:));
    filter_curr=filter{i};
    test=ismember(x_tmp,filter_curr);
    if sum(test)==length(filter_curr)
        cri(i)=1;
        index{i}=test;
    end
    
end

if sum(cri)==3
    sel.flag=1;
    sel.index=index;
else
    sel.flag=0;
    sel.index=[];
end


end

function sel=selectdata(x,y);
sel=[];
x1=x;
y1=y;
for i=1:3
    x_tmp=x(i,:);
    sel{i}=x_tmp(y{i});
    
    
end

end

