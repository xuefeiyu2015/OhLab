function Batch_OneDR(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);

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
index = 1;
for i=1:length(FilesName)
  clear OutputData;
  for j = 1:length(FilesName{i})
  
      
     load(FilesName{i}{j});
      
     
     OutputDataTemp=OutputData;
     
     Task(index)=OutputData.OneDR(1).TaskCode;
     TaskName{index}=OutputData.OneDR(1).Task;
     
     Data =OutputData.OneDR.DataStamp;

     Tgt_Good_RF(index) = Data('TgtGoodRF');
     Tgt_Bad_RF(index) = Data('TgtBad_RF');   
     p_Tgt_RF(index) = Data('p_Tgt_RF');


     PreTgtGood_RF(index) = Data('PreTgtGood_RF');
     PreTgtBad_RF(index) = Data('PreTgtBad_RF');
     p_preTgt_RF(index) = Data('p_preTgt_RF');



     

     PSTH_Good_TgtOn_RF{index} = Data('PSTH_Good_TgtOn_RF');
     PSTH_Bad_TgtOn_RF{index} = Data('PSTH_Bad_TgtOn_RF');
     PSTH_Good_TgtOn_nonRF{index} = Data('PSTH_Good_TgtOn_nonRF');
     PSTH_Bad_TgtOn_nonRF{index} = Data('PSTH_Bad_TgtOn_nonRF');

     PSTH_SacGood_nonRF{index} = Data('PSTH_SacGood_nonRF');
     PSTH_SacBad_nonRF{index} = Data('PSTH_SacBad_nonRF');
     PSTH_SacGood_RF{index} = Data('PSTH_SacGood_RF');
     PSTH_SacBad_RF{index} = Data('PSTH_SacBad_RF');

     PSTH_Time_TgtOn = Data('Time_TgtOn');
     PSTH_Time_SacOn = Data('Time_SacOn');

     %Change into prefered psth
     Diff_TgtOn_curr = Tgt_Good_RF(index) - Tgt_Bad_RF(index);
     Diff_preTgtOn_curr = PreTgtGood_RF(index) - PreTgtBad_RF(index);

     if p_Tgt_RF(index) < 0.05
         if Diff_TgtOn_curr > 0
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacGood_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacBad_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};




         else
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacBad_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacGood_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};



         end
     elseif p_preTgt_RF(index) < 0.05
         if Diff_preTgtOn_curr > 0
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacGood_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacBad_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};


         else
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacBad_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacGood_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};



         end
     else
         PSTH_Prefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};
         PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};

         PSTH_Prefer_SacOn_RF{index} = PSTH_SacGood_RF{index};
         PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacBad_RF{index};

         PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};
         PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};

         PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};
         PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};


     end


     index = index +1;



    
  end
     
     
     
end

Diff_TgtOn = p_Tgt_RF < 0.05;
Diff_PreTgtOn = p_preTgt_RF < 0.05;

SigDiff = Diff_TgtOn|Diff_PreTgtOn;




PSTH_Prefer_TgtOn_RF = cell2mat(PSTH_Prefer_TgtOn_RF');
PSTH_nonPrefer_TgtOn_RF = cell2mat(PSTH_nonPrefer_TgtOn_RF');

PSTH_Prefer_TgtOn_RF_Mean = nanmean(PSTH_Prefer_TgtOn_RF(SigDiff,:),1);
PSTH_nonPrefer_TgtOn_RF_Mean = nanmean(PSTH_nonPrefer_TgtOn_RF(SigDiff,:),1);

PSTH_Prefer_TgtOn_RF_Sem = nanstd(PSTH_Prefer_TgtOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_TgtOn_RF_Sem = nanstd(PSTH_nonPrefer_TgtOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));


PSTH_Prefer_SacOn_RF = cell2mat(PSTH_Prefer_SacOn_RF');
PSTH_nonPrefer_SacOn_RF = cell2mat(PSTH_nonPrefer_SacOn_RF');

PSTH_Prefer_SacOn_RF_Mean = nanmean(PSTH_Prefer_SacOn_RF(SigDiff,:),1);
PSTH_nonPrefer_SacOn_RF_Mean = nanmean(PSTH_nonPrefer_SacOn_RF(SigDiff,:),1);

PSTH_Prefer_SacOn_RF_Sem = nanstd(PSTH_Prefer_SacOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_SacOn_RF_Sem = nanstd(PSTH_nonPrefer_SacOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));


PSTH_Prefer_TgtOn_nonRF = cell2mat(PSTH_Prefer_TgtOn_nonRF');
PSTH_nonPrefer_TgtOn_nonRF = cell2mat(PSTH_nonPrefer_TgtOn_nonRF');

PSTH_Prefer_TgtOn_nonRF_Mean = nanmean(PSTH_Prefer_TgtOn_nonRF(SigDiff,:),1);
PSTH_nonPrefer_TgtOn_nonRF_Mean = nanmean(PSTH_nonPrefer_TgtOn_nonRF(SigDiff,:),1);

PSTH_Prefer_TgtOn_nonRF_Sem = nanstd(PSTH_Prefer_TgtOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_TgtOn_nonRF_Sem = nanstd(PSTH_nonPrefer_TgtOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));


PSTH_Prefer_SacOn_nonRF = cell2mat(PSTH_Prefer_SacOn_nonRF');
PSTH_nonPrefer_SacOn_nonRF = cell2mat(PSTH_nonPrefer_SacOn_nonRF');

PSTH_Prefer_SacOn_nonRF_Mean = nanmean(PSTH_Prefer_SacOn_nonRF(SigDiff,:),1);
PSTH_nonPrefer_SacOn_nonRF_Mean = nanmean(PSTH_nonPrefer_SacOn_nonRF(SigDiff,:),1);

PSTH_Prefer_SacOn_nonRF_Sem = nanstd(PSTH_Prefer_SacOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_SacOn_nonRF_Sem = nanstd(PSTH_nonPrefer_SacOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));




if ShowFigureFlag
    
FigureStartNum=100;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_PrevsNonpre';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
subplot(2,2,1)

plot(PSTH_Time_TgtOn,PSTH_Prefer_TgtOn_RF_Mean,'-r');
hold on
plot(PSTH_Time_TgtOn,PSTH_nonPrefer_TgtOn_RF_Mean,'-b');

title('Contralateral(TgtOn)')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,160]);
box off

 subplot(2,2,2)


 plot(PSTH_Time_SacOn,PSTH_Prefer_SacOn_RF_Mean,'-r');
hold on
plot(PSTH_Time_SacOn,PSTH_nonPrefer_SacOn_RF_Mean,'-b');
title('Contralateral(SacOn)')

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylim([-2,40])
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

subplot(2,2,3)

plot(PSTH_Time_TgtOn,PSTH_Prefer_TgtOn_nonRF_Mean,'-r');
hold on
plot(PSTH_Time_TgtOn,PSTH_nonPrefer_TgtOn_nonRF_Mean,'-b');

title('Ipsilateral(TgtOn)')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,160]);
box off

 subplot(2,2,4)


 plot(PSTH_Time_SacOn,PSTH_Prefer_SacOn_nonRF_Mean,'-r');
hold on
plot(PSTH_Time_SacOn,PSTH_nonPrefer_SacOn_nonRF_Mean,'-b');
title('Ipsilateral(SacOn)')

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylim([-2,40])
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

%}
end







end