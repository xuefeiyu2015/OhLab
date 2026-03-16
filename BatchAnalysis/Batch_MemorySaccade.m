function Batch_MemorySaccade(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
%{
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
%}
%Load file
for i=1:length(FilesName)
  clear OutputData;
  
     load(FilesName{i}{1});
     
     OutputDataTemp=OutputData;
  %   Task(i)=OutputData.MemorySaccade(1).TaskCode;
  %   TaskName{i}=OutputData.MemorySaccade(1).Task;
     
  %   Data =OutputData.MemorySaccade(1).DataStamp;

   Task(i)=OutputData.MemorySaccadeLatency(1).TaskCode;
     TaskName{i}=OutputData.MemorySaccadeLatency(1).Task;
     
     Data =OutputData.MemorySaccadeLatency(1).DataStamp;

     PSTH_Tgt_RF{i} = Data('PSTH_TgtOn_RF');
     PSTH_Sac_RF{i} = Data('PSTH_SacOn_RF');

   %  RF_h_v(i) = Data('RF_h_v');
   %  RF_h_m(i) = Data('RF_h_m');
%
     CellType(i) = Data('CellType');




    
     
     
end
PSTH_Time_TgtOn = Data('PSTH_Time_TgtOn');
PSTH_Time_SacOn = Data('PSTH_Time_SacOn');



PSTH_Tgt_RF = cell2mat(cellfun(@(x) x(40:150), PSTH_Tgt_RF,'UniformOutput',false)');
PSTH_Sac_RF = cell2mat(PSTH_Sac_RF');

PSTH_Time_TgtOn = PSTH_Time_TgtOn(40:150);


CellTypeUnique = unique(CellType);
for i = 1:length(CellTypeUnique)
    ct = CellType == CellTypeUnique(i);
    


        PSTH_Tgt_Mean(i,:) = mean(PSTH_Tgt_RF(ct,:),1);
        PSTH_Sac_Mean(i,:) = mean(PSTH_Sac_RF(ct,:),1);

        PSTH_Tgt_Sem(i,:) = std(PSTH_Tgt_RF(ct,:),[],1)/sqrt(sum(ct));
        PSTH_Sac_Sem(i,:) = std(PSTH_Sac_RF(ct,:),[],1)/sqrt(sum(ct));




end

%Remove none type
UCType = CellTypeUnique(CellTypeUnique > 0);


PSTH_Tgt_Mean = PSTH_Tgt_Mean(CellTypeUnique > 0,:);
PSTH_Sac_Mean = PSTH_Sac_Mean(CellTypeUnique > 0,:);

if ShowFigureFlag
 
FigureStartNum=100;
FigureIndex=1;

figtitlestr{FigureIndex}='PSTH_MemorySac_Type';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

for i = 1:length(UCType)
   subplot(length(UCType),2,(i-1)*2+1)
   plot(PSTH_Time_TgtOn,PSTH_Tgt_Mean(i,:),'-r');
   shadedErrorBar(PSTH_Time_TgtOn,PSTH_Tgt_Mean(i,:),PSTH_Tgt_Sem(i,:),'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
 xlabel('Time from target onset(ms)');
 box off;
 % ylim([0,100]);
 

   subplot(length(UCType),2,(i-1)*2+2)
    shadedErrorBar(PSTH_Time_SacOn,PSTH_Sac_Mean(i,:),PSTH_Sac_Sem(i,:),'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
   xlabel('Time from saccade onset(ms)');
   box off;
 %ylim([0,100]);
 xlim([-500,200]);


end



end
    