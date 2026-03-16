function Batch_OptoOneDR(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);

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

     if isfield(OutputData,'OneDR')
     
     file_use{index} = FilesName{i}{j};
     
     Data =OutputData.FreeViewNeural.DataStamp;

     PSTH_Stim_Mean{index} = Data('PSTH_Stim_Mean');
     
     AS(index) = Data('AcStrength');
     Latency(index) = Data('ResponseLatency');
%     p_AS(index) = Data('Ac_p');
     PMS(index) = Data('PostModulationStrength');
     PMS_p(index) = Data('PostModulation_p');
     SM(index) = Data('StimModulation');
     SM_p(index) = Data('StimModulation_p');
     PSTH_First{index} = Data('PSTH_FirstHalf');
     PSTH_Second{index} = Data('PSTH_SecondHalf');
     AS_First(index) = Data('AS_First');
     AS_Second(index) = Data('AS_Second');
    PSTH_Time = Data('PSTH_Time');
          
      
     
   
     
     Data =OutputData.OneDR.DataStamp;

     Tgt_Good_RF(index) = Data('TgtGoodRF');
     Tgt_Bad_RF(index) = Data('TgtBad_RF');   
     p_Tgt_RF(index) = Data('p_Tgt_RF');


     PreTgtGood_RF(index) = Data('PreTgtGood_RF');
     PreTgtBad_RF(index) = Data('PreTgtBad_RF');
     p_preTgt_RF(index) = Data('p_preTgt_RF');



     
%{
     PSTH_Good_TgtOn_RF{index} = Data('PSTH_Good_TgtOn_RF');
     PSTH_Bad_TgtOn_RF{index} = Data('PSTH_Bad_TgtOn_RF');
     PSTH_Good_TgtOn_nonRF{index} = Data('PSTH_Good_TgtOn_nonRF');
     PSTH_Bad_TgtOn_nonRF{index} = Data('PSTH_Bad_TgtOn_nonRF');

     PSTH_SacGood_nonRF{index} = Data('PSTH_SacGood_nonRF');
     PSTH_SacBad_nonRF{index} = Data('PSTH_SacBad_nonRF');
     PSTH_SacGood_RF{index} = Data('PSTH_SacGood_RF');
     PSTH_SacBad_RF{index} = Data('PSTH_SacBad_RF');
%}
      PSTH_Good_TgtOn_RF{index} = Data('PSTH_Good_TgtOn_RF_z');
     PSTH_Bad_TgtOn_RF{index} = Data('PSTH_Bad_TgtOn_RF_z');
     PSTH_Good_TgtOn_nonRF{index} = Data('PSTH_Good_TgtOn_nonRF_z');
     PSTH_Bad_TgtOn_nonRF{index} = Data('PSTH_Bad_TgtOn_nonRF_z');

     PSTH_SacGood_nonRF{index} = Data('PSTH_SacGood_nonRF_z');
     PSTH_SacBad_nonRF{index} = Data('PSTH_SacBad_nonRF_z');
     PSTH_SacGood_RF{index} = Data('PSTH_SacGood_RF_z');
     PSTH_SacBad_RF{index} = Data('PSTH_SacBad_RF_z');

%{
'PSTH_Good_TgtOn_RF_z','PSTH_Bad_TgtOn_RF_z','PSTH_Good_TgtOn_nonRF_z','PSTH_Bad_TgtOn_nonRF_z',...
    'PSTH_SacGood_RF_z','PSTH_SacBad_RF_z',...
    'TgtGood_RF_z','TgtBad_RF_z','PreTgtGood_RF_z','PreTgtBad_RF_z'
%}

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


          
     end
     index = index +1;
     end % Endof if isfield 


    

  end
     
     
     
end

Diff_TgtOn = p_Tgt_RF < 0.05;
Diff_PreTgtOn = p_preTgt_RF < 0.05;

%SigDiff = (Diff_TgtOn|Diff_PreTgtOn)&SM_p<0.05;
%SigDiff = SM_p>0.05;
Activation_index = SM>0;
Inhibition_index = SM<=0;

SigDiff = SM_p<0.05 & Latency<=20 ;
NonSigDiff = SM_p>=0.05;

N_tag = sum(SigDiff)
N_nontag = sum(NonSigDiff)




PSTH_Prefer_TgtOn_RF = cell2mat(PSTH_Prefer_TgtOn_RF');
PSTH_nonPrefer_TgtOn_RF = cell2mat(PSTH_nonPrefer_TgtOn_RF');

PSTH_Prefer_SacOn_nonRF = cell2mat(PSTH_Prefer_SacOn_nonRF');
PSTH_nonPrefer_SacOn_nonRF = cell2mat(PSTH_nonPrefer_SacOn_nonRF');

PSTH_Prefer_SacOn_RF = cell2mat(PSTH_Prefer_SacOn_RF');
PSTH_nonPrefer_SacOn_RF = cell2mat(PSTH_nonPrefer_SacOn_RF');

PSTH_Prefer_TgtOn_nonRF = cell2mat(PSTH_Prefer_TgtOn_nonRF');
PSTH_nonPrefer_TgtOn_nonRF = cell2mat(PSTH_nonPrefer_TgtOn_nonRF');

%For opto-tagged neurons




PSTH_Prefer_TgtOn_RF_Mean = nanmean(PSTH_Prefer_TgtOn_RF(SigDiff,:),1);
PSTH_nonPrefer_TgtOn_RF_Mean = nanmean(PSTH_nonPrefer_TgtOn_RF(SigDiff,:),1);

PSTH_Prefer_TgtOn_RF_Sem = nanstd(PSTH_Prefer_TgtOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_TgtOn_RF_Sem = nanstd(PSTH_nonPrefer_TgtOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));






PSTH_Prefer_SacOn_RF_Mean = nanmean(PSTH_Prefer_SacOn_RF(SigDiff,:),1);
PSTH_nonPrefer_SacOn_RF_Mean = nanmean(PSTH_nonPrefer_SacOn_RF(SigDiff,:),1);

PSTH_Prefer_SacOn_RF_Sem = nanstd(PSTH_Prefer_SacOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_SacOn_RF_Sem = nanstd(PSTH_nonPrefer_SacOn_RF(SigDiff,:),[],1)/sqrt(sum(SigDiff));





PSTH_Prefer_TgtOn_nonRF_Mean = nanmean(PSTH_Prefer_TgtOn_nonRF(SigDiff,:),1);
PSTH_nonPrefer_TgtOn_nonRF_Mean = nanmean(PSTH_nonPrefer_TgtOn_nonRF(SigDiff,:),1);

PSTH_Prefer_TgtOn_nonRF_Sem = nanstd(PSTH_Prefer_TgtOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_TgtOn_nonRF_Sem = nanstd(PSTH_nonPrefer_TgtOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));





PSTH_Prefer_SacOn_nonRF_Mean = nanmean(PSTH_Prefer_SacOn_nonRF(SigDiff,:),1);
PSTH_nonPrefer_SacOn_nonRF_Mean = nanmean(PSTH_nonPrefer_SacOn_nonRF(SigDiff,:),1);

PSTH_Prefer_SacOn_nonRF_Sem = nanstd(PSTH_Prefer_SacOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));
PSTH_nonPrefer_SacOn_nonRF_Sem = nanstd(PSTH_nonPrefer_SacOn_nonRF(SigDiff,:),[],1)/sqrt(sum(SigDiff));




%For nonopto-tagged neurons

PSTH_Prefer_TgtOn_RF_Mean_nontag = nanmean(PSTH_Prefer_TgtOn_RF(NonSigDiff,:),1);
PSTH_nonPrefer_TgtOn_RF_Mean_nontag = nanmean(PSTH_nonPrefer_TgtOn_RF(NonSigDiff,:),1);

PSTH_Prefer_TgtOn_RF_Sem_nontag = nanstd(PSTH_Prefer_TgtOn_RF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));
PSTH_nonPrefer_TgtOn_RF_Sem_nontag = nanstd(PSTH_nonPrefer_TgtOn_RF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));

%PSTH_Prefer_TgtOn_RF_Sem_nontag = nanstd(PSTH_Prefer_TgtOn_RF(NonSigDiff,:),[],1)/(sum(NonSigDiff));
%PSTH_nonPrefer_TgtOn_RF_Sem_nontag = nanstd(PSTH_nonPrefer_TgtOn_RF(NonSigDiff,:),[],1)/(sum(NonSigDiff));




PSTH_Prefer_SacOn_RF_Mean_nontag = nanmean(PSTH_Prefer_SacOn_RF(NonSigDiff,:),1);
PSTH_nonPrefer_SacOn_RF_Mean_nontag = nanmean(PSTH_nonPrefer_SacOn_RF(NonSigDiff,:),1);

PSTH_Prefer_SacOn_RF_Sem_nontag = nanstd(PSTH_Prefer_SacOn_RF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));
PSTH_nonPrefer_SacOn_RF_Sem_nontag = nanstd(PSTH_nonPrefer_SacOn_RF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));

%PSTH_Prefer_SacOn_RF_Sem_nontag = nanstd(PSTH_Prefer_SacOn_RF(NonSigDiff,:),[],1)/(sum(NonSigDiff));
%PSTH_nonPrefer_SacOn_RF_Sem_nontag = nanstd(PSTH_nonPrefer_SacOn_RF(NonSigDiff,:),[],1)/(sum(NonSigDiff));




PSTH_Prefer_TgtOn_nonRF_Mean_nontag = nanmean(PSTH_Prefer_TgtOn_nonRF(NonSigDiff,:),1);
PSTH_nonPrefer_TgtOn_nonRF_Mean_nontag = nanmean(PSTH_nonPrefer_TgtOn_nonRF(NonSigDiff,:),1);

PSTH_Prefer_TgtOn_nonRF_Sem_nontag = nanstd(PSTH_Prefer_TgtOn_nonRF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));
PSTH_nonPrefer_TgtOn_nonRF_Sem_nontag = nanstd(PSTH_nonPrefer_TgtOn_nonRF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));

%PSTH_Prefer_TgtOn_nonRF_Sem_nontag = nanstd(PSTH_Prefer_TgtOn_nonRF(NonSigDiff,:),[],1)/(sum(NonSigDiff));
%PSTH_nonPrefer_TgtOn_nonRF_Sem_nontag = nanstd(PSTH_nonPrefer_TgtOn_nonRF(NonSigDiff,:),[],1)/(sum(NonSigDiff));




PSTH_Prefer_SacOn_nonRF_Mean_nontag = nanmean(PSTH_Prefer_SacOn_nonRF(NonSigDiff,:),1);
PSTH_nonPrefer_SacOn_nonRF_Mean_nontag = nanmean(PSTH_nonPrefer_SacOn_nonRF(NonSigDiff,:),1);

PSTH_Prefer_SacOn_nonRF_Sem_nontag = nanstd(PSTH_Prefer_SacOn_nonRF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));
PSTH_nonPrefer_SacOn_nonRF_Sem_nontag = nanstd(PSTH_nonPrefer_SacOn_nonRF(NonSigDiff,:),[],1)/sqrt(sum(NonSigDiff));
%PSTH_Prefer_SacOn_nonRF_Sem_nontag = nanstd(PSTH_Prefer_SacOn_nonRF(NonSigDiff,:),[],1)/(sum(NonSigDiff));
%PSTH_nonPrefer_SacOn_nonRF_Sem_nontag = nanstd(PSTH_nonPrefer_SacOn_nonRF(NonSigDiff,:),[],1)/(sum(NonSigDiff));

%Distribution
SigLoc_Taged = (Diff_TgtOn|Diff_PreTgtOn)&SigDiff ;
NonSigLoc_Taged = (~(Diff_TgtOn|Diff_PreTgtOn))&SigDiff ;
Prop_Sig_Taged = sum(SigLoc_Taged)/sum(SigDiff );
Prop_NonSig_Taged = sum(NonSigLoc_Taged)/sum(SigDiff );

Num_Sig_Taged = sum(SigLoc_Taged);
Num_NonSig_Taged = sum(NonSigLoc_Taged);


SigLoc_NonTaged = (Diff_TgtOn|Diff_PreTgtOn)&NonSigDiff ;
NonSigLoc_NonTaged = (~(Diff_TgtOn|Diff_PreTgtOn))&NonSigDiff ;
Prop_Sig_NonTaged = sum(SigLoc_NonTaged)/sum(NonSigDiff );
Prop_NonSig_NonTaged = sum(NonSigLoc_NonTaged)/sum(NonSigDiff );

Num_Sig_NonTaged = sum(SigLoc_NonTaged);
Num_NonSig_NonTaged = sum(NonSigLoc_NonTaged);


[chi_squared, df, p_prop, is_significant] = chisq_prop_test(Num_Sig_Taged , Num_NonSig_Taged, Num_Sig_NonTaged , Num_NonSig_NonTaged);
disp(sprintf('p value for the difference between the two proportion is :%1.2f',p_prop));


%Mat_dist = [Prop_Sig_Taged,Prop_Sig_NonTaged;Prop_NonSig_Taged,Prop_NonSig_NonTaged];
Mat_dist = [Num_Sig_Taged,Num_Sig_NonTaged;Num_NonSig_Taged,Num_NonSig_NonTaged];



ActivationStrength = SM(SigLoc_Taged);
Latency_sig =Latency(SigLoc_Taged);


Pro_E_tag = sum(ActivationStrength > 0)/length(ActivationStrength);
Pro_I_tag = sum(ActivationStrength <0)/length(ActivationStrength);




Latency_nonsig =Latency(NonSigLoc_Taged);








if ShowFigureFlag
    
FigureStartNum=100;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_PrevsNonpre_OptoTaged';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
ax1 = subplot(2,2,1);

shadedErrorBar(PSTH_Time_TgtOn,PSTH_Prefer_TgtOn_RF_Mean,PSTH_Prefer_TgtOn_RF_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_TgtOn,PSTH_nonPrefer_TgtOn_RF_Mean,PSTH_nonPrefer_TgtOn_RF_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

title('Contralateral(TgtOn)')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,160]);
box off

 ax2 = subplot(2,2,2);


 shadedErrorBar(PSTH_Time_SacOn,PSTH_Prefer_SacOn_RF_Mean,PSTH_Prefer_SacOn_RF_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_SacOn,PSTH_nonPrefer_SacOn_RF_Mean,PSTH_nonPrefer_SacOn_RF_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);
title('Contralateral(SacOn)')

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-2,40])
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

ax3 = subplot(2,2,3);

shadedErrorBar(PSTH_Time_TgtOn,PSTH_Prefer_TgtOn_nonRF_Mean,PSTH_Prefer_TgtOn_nonRF_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_TgtOn,PSTH_nonPrefer_TgtOn_nonRF_Mean,PSTH_nonPrefer_TgtOn_nonRF_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

title('Ipsilateral(TgtOn)')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,160]);
box off

 ax4 = subplot(2,2,4);


shadedErrorBar(PSTH_Time_SacOn,PSTH_Prefer_SacOn_nonRF_Mean,PSTH_Prefer_SacOn_nonRF_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_SacOn,PSTH_nonPrefer_SacOn_nonRF_Mean,PSTH_nonPrefer_SacOn_nonRF_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);
title('Ipsilateral(SacOn)')

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-2,40])
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

linkaxes([ax1 ax2 ax3 ax4],'y');

%}
%%
FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;
    
figtitlestr{FigureIndex}='PSTH_PrevsNonpre_NotOptoTaged';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
ax1 = subplot(2,2,1);

shadedErrorBar(PSTH_Time_TgtOn,PSTH_Prefer_TgtOn_RF_Mean_nontag,PSTH_Prefer_TgtOn_RF_Sem_nontag,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_TgtOn,PSTH_nonPrefer_TgtOn_RF_Mean_nontag,PSTH_nonPrefer_TgtOn_RF_Sem_nontag,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

title('Contralateral(TgtOn)')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,160]);
box off

 ax2 = subplot(2,2,2);


 shadedErrorBar(PSTH_Time_SacOn,PSTH_Prefer_SacOn_RF_Mean_nontag,PSTH_Prefer_SacOn_RF_Sem_nontag,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_SacOn,PSTH_nonPrefer_SacOn_RF_Mean_nontag,PSTH_nonPrefer_SacOn_RF_Sem_nontag,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);
title('Contralateral(SacOn)')

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-2,40])
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

ax3 = subplot(2,2,3);

shadedErrorBar(PSTH_Time_TgtOn,PSTH_Prefer_TgtOn_nonRF_Mean_nontag,PSTH_Prefer_TgtOn_nonRF_Sem_nontag,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_TgtOn,PSTH_nonPrefer_TgtOn_nonRF_Mean_nontag,PSTH_nonPrefer_TgtOn_nonRF_Sem_nontag,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

title('Ipsilateral(TgtOn)')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,160]);
box off

 ax4 = subplot(2,2,4);


 shadedErrorBar(PSTH_Time_SacOn,PSTH_Prefer_SacOn_nonRF_Mean_nontag,PSTH_Prefer_SacOn_nonRF_Sem_nontag,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_SacOn,PSTH_nonPrefer_SacOn_nonRF_Mean_nontag,PSTH_nonPrefer_SacOn_nonRF_Sem_nontag,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);
title('Ipsilateral(SacOn)')

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-2,40])
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

linkaxes([ax1 ax2 ax3 ax4],'y');

%}
%% 
FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;
    
figtitlestr{FigureIndex}='Summary';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
subplot(2,2,1);
b = bar(Mat_dist','stacked');
b(1).FaceColor = 'r';
b(2).FaceColor = 'b';
xticks([1,2]);
xticklabels({'Tag','Non-tag'});

legend({'Sig','Non-sig'});
ylabel('Number of neurons');

box off
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');


subplot(2,2,2);
b = bar([1,2],[Pro_E_tag,Pro_I_tag]);
b.FaceColor = 'flat';
b.CData =[[1,0,0];[0,0,1]]; 

xticks([1,2]);
xticklabels({'Exc','Inh'});
yticks([0,0.5,1]);
ylabel('Proportion');

box off
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

subplot(2,2,3)
histogram(Latency_sig,'BinWidth',10,'FaceColor','r','LineWidth',3);
xlim([-10,100]);
ylabel('Number of neurons');
xlabel('Latency of Opto-stim');
title('Sig Loc Value');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off;

subplot(2,2,4)
histogram(Latency_nonsig,'BinWidth',10,'FaceColor','k','LineWidth',3);
ylabel('Number of neurons');
xlim([-10,100]);
ylabel('Number of neurons');
xlabel('Latency of Opto-stim');
title('NonSig Loc Value');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off;

end



end