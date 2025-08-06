
function Batch_Optogenetics(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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


FilesName_Each=Data.FileName;
ChannelID=Data.ChannelNumber;

index =1;

ChannelFull=[];

FileIndex=[];
AbsoluteIndex = [];


for i=1:length(FilesName)
  clear OutputData;
   ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];

    for j = 1:length(FilesName{i})

        %For Video Free View
  
     load(FilesName{i}{j});
     
     OutputDataTemp=OutputData;

     FilesNameAll{index}=FilesName_Each{i};

     AbsoluteIndex(index)=index;

     if index == 1
         FileIndex(index) = 1;

     else
         if ~strcmp(FilesNameAll{index},FilesNameAll{index-1})
             FileIndex(index) = FileIndex(index-1)+1;
         else
             FileIndex(index) = FileIndex(index-1);

         end
     end





      Data =OutputData.FreeViewNeural.DataStamp;
%
try

    PSTH_Stim_Norm{index} = Data('PSTH_Stim_Mean_norm');
catch
    keyboard
end
     PSTH_Control_Norm{index} = Data('PSTH_Control_Mean_norm');
     PSTH_Stim_Mean_z{index} = Data('PSTH_Stim_Mean_z');

   %   PSTH_Stim_Mean_z{index} = Data('PSTH_Stim_Mean_z');%Use the z scored ones
   %  PSTH_Control_Mean_z{index} = Data('PSTH_Control_Mean_z');%Use the z scored ones
%}
  %    PSTH_Stim_Norm{index} = Data('PSTH_Stim_Mean_z');
  %   PSTH_Control_Norm{index} = Data('PSTH_Control_Mean_z');

      PSTH_Stim_Mean{index} = Data('PSTH_Stim_Mean');
     PSTH_Control_Mean{index} = Data('PSTH_Control_Mean');

     FR_Mean_Video(index) = Data('FR_Mean');

     p_StimControl(index) = Data('p_StimControl');

     p_Cri20(index) = Data('p_cri20');
     p_Cri50(index) = Data('p_cri50');
     Latency2(index) = Data('Latency2');

     PSTH_Stim_Norm2{index} = Data('psth_stim_normed2');


   


     
     AS(index) = Data('AcDur');
     Latency(index) = Data('ResponseLatency');
     p_AS(index) = Data('Ac_p');
     PMS(index) = Data('PostModulationStrength');
     PMS_p(index) = Data('PostModulation_p');
     SM(index) = Data('StimModulation');
     SM_p(index) = Data('StimModulation_p');
     %PSTH_First{index} = Data('PSTH_FirstHalf');
     %PSTH_Second{index} = Data('PSTH_SecondHalf');
     PSTH_First{index} = Data('PSTH_Stim_Mean_First_z');
     PSTH_Second{index} = Data('PSTH_Stim_Mean_Second_z');

     p_whole_cri50{index} = Data('p_whole_cri50');



     


     AS_First(index) = Data('AS_First');
     AS_Second(index) = Data('AS_Second');

     PSTH_Opto_Time = Data('PSTH_Time');

 
     index = index+1;   
    end %End of for file name
     
end


%Plot controller
PlotOptoEffect = 1;


%Analysis controller 
AnalysisOS = 1;% For Video Free View


%Selection 
%ManuelOut
CellIndex = 1:length(FR_Mean_Video);
%ManuelCheck=[550 469	478	477	473	479	425	408	353	351	320	239	240	217	213	201	178	163	160	139	141	68	43	48 115 132 137 157 181 200 221 496 516 529 538 581 623 97 131 159 315 325 329 343 350 400 427 435];
ManuelCheck=[496];
ManuelOut = ~ismember(CellIndex,ManuelCheck);


Cri1 = FR_Mean_Video>1;



ScreenedIndex = CellIndex(Cri1 & ManuelOut);

ScreenCell = ismember(CellIndex,ScreenedIndex);

ScreenCellNo = CellIndex(ScreenCell);


%% For VideoFreeView

PSTH = cell2mat(PSTH_Stim_Mean(ScreenCell)'); 
PSTH_control = cell2mat(PSTH_Control_Mean(ScreenCell)');

PSTH_norm = cell2mat(PSTH_Stim_Norm(ScreenCell)'); 
PSTH_control_norm = cell2mat(PSTH_Control_Norm(ScreenCell)');


PSTH_norm2 =cell2mat(PSTH_Stim_Norm2(ScreenCell)');



PSTH_z = cell2mat(PSTH_Stim_Mean_z(ScreenCell)');


p_whole_cri50 = cell2mat(p_whole_cri50(ScreenCell));

Latency=Latency(ScreenCell);

SM_p=SM_p(ScreenCell);
SM=SM(ScreenCell);
PMS_p=PMS_p(ScreenCell);

%maxnum_ms = max(cellfun(@numel,PSTH_Tgt_RF));







%{
 p_Tgt_RF_stable(index) = Data('p_Tgt_RF');


      %  PreTgtGood_RF(index) = Data('PreTgtGood_RF');
      %  PreTgtBad_RF(index) = Data('PreTgtBad_RF');
        p_Sac_RF_stable(index) = Data('p_Sac_RF');

  TgtGood_RF_Stable(index) = NaN;
        TgtBad_RF_Stable(index) = NaN;

         TgtGood_nonRF_Stable(index) = NaN;
        TgtBad_nonRF_Stable(index) = NaN;

        SacGood_RF_Stable(index) = NaN;
        SacBad_RF_Stable(index) = NaN;

         SacGood_nonRF_Stable(index) = NaN;
        SacBad_nonRF_Stable(index) = NaN;

        Value_Vis_RF(index) = NaN;

        Value_Vis_nonRF(index) = NaN;

        Value_Sac_RF(index) = NaN;

        Value_Sac_nonRF(index) = NaN;


%}

%PSTH_Tgt_RF_mat = cellfun(@(x) [x,NaN*ones(1,maxnum_ms-length(x))],PSTH_Tgt_RF,'UniformOutput',0);

%%
%Separate cells into different modulation groups according to laser
%response


%Separate into 5 groups 


Fast_M = Latency<5;
Slow_M = Latency>=5;
Excitation = SM> 0;
Inhibition = SM<=0;

%Sig_m = SM_p < 0.05;
%Sig_m = ~isnan(Latency) & AS >=5;

Sig_m = (p_whole_cri50<0.05 | SM_p<0.05 | PMS_p<0.05)&(~isnan(Latency));
%Sig_m = p_whole_cri50>-0.1;
%Latency
Non_Sig_m = ~Sig_m;

SepLevel = 0.5;

Non_Sig_level2 = p_whole_cri50>SepLevel & SM_p>SepLevel;
Non_Sig_level1 = ~Non_Sig_level2  & Non_Sig_m;

%Non_Sig_m = Non_Sig_level2;


Ex_Fast_Sig = Fast_M & Excitation & Sig_m;

Ex_Slow_Sig = Slow_M & Excitation & Sig_m;


%In_Fast_Sig = Fast_M & Inhibition & Sig_m;

%In_Slow_Sig = Slow_M & Inhibition & Sig_m;

In_Sig = Inhibition & Sig_m;

%PSTH_Use = PSTH_norm;
%PSTH_Use = PSTH; 
%PSTH_Use = PSTH_z; 
PSTH_Use = PSTH_norm2;

ModuType = NaN*ones(1,length(Sig_m));
NSig_Level = NaN*ones(1,length(Sig_m));



%%
ModuType(Ex_Fast_Sig) =1;
ModuType(Ex_Slow_Sig) =2;

ModuType(In_Sig) =3;
%ModuType(In_Slow_Sig) =4;

ModuType(Non_Sig_m) =4;

NSig_Level(Non_Sig_level1&Excitation) =1;
NSig_Level(Non_Sig_level1&Inhibition) =2;

NSig_Level(Non_Sig_level2&Excitation) =3;
NSig_Level(Non_Sig_level2&Inhibition) =4;

NSig_Level=NSig_Level';


ModulString = {'Fast Excitation','Slow Excitation','Inhibition','No Effect'};

NonSigString = {'NonSig_Ex_l1','NonSig_In_l1','NonSig_Ex_l2','NonSig_In_l2'};

Latency_Ex = Latency(Ex_Fast_Sig|Ex_Slow_Sig);
Latency_In = Latency(In_Sig);



UniqueMol = unique(ModuType);

UniqueMol = UniqueMol(~isnan(UniqueMol));

ModuType=ModuType';%Transfer to column;

UniqueNonSig = unique(NSig_Level);
UniqueNonSig=UniqueNonSig(~isnan(UniqueNonSig));

%For pie plot
ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);
FastExcitationNum = sum(Ex_Fast_Sig);
SlowExcitationNum = sum(Ex_Slow_Sig);



%

%%
%{
%Check = isnan(Latency) & Sig_m;

%A=PSTH_Use(Check,:);
A=PSTH_Use;
figure
%index = 1;
for index = 0:size(A,1)-1

    
    
    %if ModuType(i+1) == 2
    fig_index = floor(index+1/16)+1;
    subplot_index = mod(index+1,16)+1;
    figure(fig_index);
    subplot(4,4,subplot_index);

    

        plot(A(index+1,:));
      %  index=index+1;
   % end
   
end
%}
%{
figure
subplot(2,1,1)
edge = 0:2.5:150;
histogram(Latency(Excitation & Sig_m),edge );
subplot(2,1,2)
histogram(Latency(Inhibition & Sig_m),edge )
%}


%% Manuel Check for All cells 
%For each group
%{
Group1_name = AbsoluteIndex(ModuType == 1);

Group1_PSTH_Stim = PSTH_Use(ModuType == 1,:);

Group2_name = AbsoluteIndex(ModuType == 2);

Group2_PSTH_Stim = PSTH_Use(ModuType == 2,:);

Group3_name = AbsoluteIndex(ModuType == 3);

Group3_PSTH_Stim = PSTH_Use(ModuType == 3,:);

Group4_name = AbsoluteIndex(ModuType == 4);

Group4_PSTH_Stim = PSTH_Use(ModuType == 4,:);
%}

%index = 1;

%{
for i = 0:length(Group1_name )-1

    fig_index = floor(i/16)+1;
    subplot_index = mod(i,16)+1;
    figure(fig_index);
    subplot(4,4,subplot_index);
    

        plot(Group1_PSTH_Stim(i+1,:));
       % index=index+1;
       title(Group1_name(i+1))

end
%}
%{
for i = 0:length(Group2_name )-1

    fig_index = floor(i/16)+1;
    subplot_index = mod(i,16)+1;
    figure(fig_index);
    subplot(4,4,subplot_index);
    

        plot(Group2_PSTH_Stim(i+1,:));
       % index=index+1;
       title(Group2_name(i+1))

end

%}

%checkGroup(Group4_name,Group4_PSTH_Stim);

  
%AbsoluteIndex

%keyboard


for i = 1:length(UniqueMol)
    sel = ModuType == UniqueMol(i);

    ModuTypeCount(i) = sum(sel);

    PSTH_Stim_Group = PSTH_Use(sel,:);

    names = AbsoluteIndex(sel);

    
  

    PSTH_Stim_Mean_Group(i,:) = mean(PSTH_Stim_Group,1);

    PSTH_Stim_Sem_Group(i,:) = std(PSTH_Stim_Group,[],1)/sqrt(sum(sel));
    

    Latency_Mean_Group(i) = mean(Latency(sel));

    Latency_Sem_Group(i) = std(Latency(sel))/sqrt(sum(sel));

    MS_Mean_Group(i) = mean(abs(SM(sel)));
    MS_Sem_Group(i) = std(abs(SM(sel)))/sqrt(sum(sel));


end
%% If ananlyze optostimulation effect 
if AnalysisOS

%Anova for latency test for all groups 
%[p,tbl,stats] = anova1(hogg);
[p_Anova_Latency,ANOVATAB,stats]= anova1(Latency(ModuType<4),ModuType(ModuType<4),'off');

disp(sprintf('Anova test for latency for groups:%d',p_Anova_Latency));

%Perform multiple comparison
results = multcompare(stats);

disp('Mutiple comparison between different groups for latency: ');
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

%Anova for modulation strength
[p_Anova_Strength,ANOVATAB,stats]= anova1(abs(SM(ModuType<4)),ModuType(ModuType<4),'off');

disp(sprintf('Anova test for strength for groups:%d',p_Anova_Strength));

%Perform multiple comparison
results = multcompare(stats);

disp('Mutiple comparison between different groups for strength: ');
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

close all;
end 
%}

%% For nonsig neurons 
for i = 1:length(UniqueNonSig)
    sel = NSig_Level == UniqueNonSig(i);
    NonSigCount(i) = sum(sel);


    PSTH_NonSigStim_Group = PSTH_Use(sel,:);

   

  

    PSTH_NonSigStim_Mean_Group(i,:) = mean(PSTH_NonSigStim_Group,1);

    PSTH_NonSigStim_Sem_Group(i,:) = std(PSTH_NonSigStim_Group,[],1)/sqrt(sum(sel));

   



end



%{
%Onebyone check
for i = 0:size(PSTH_Use,1)-1

    fig_index = floor(i/16)+1;
    subplot_index = mod(i,16)+1;


    figure(fig_index);

    subplot(4,4,subplot_index);
    if ModuType(i+1) == 4
    plot(PSTH_Use(i+1,:));
    title(sprintf('%d,Type:%d',ScreenCellNo(i+1),ModuType(i+1)));
    end

end
%}





    
%% One by One check

%%Check the ROC
%{
 for i = 1:length(UniqueMol)
    % AUC_Stable_RF_Type{i}=AUC_Stable_RF(sel);
     %   AUC_Stable_nonRF_Type{i}=AUC_Stable_nonRF(sel);
    p_RF =  p_Tgt_RF_stable_Mo{i};
    p_nonRF = p_Tgt_nonRF_stable_Mo{i};

    auc_RF = AUC_Stable_RF_Type{i}';
    auc_nonRF = AUC_Stable_nonRF_Type{i}';

    sel_good_prefer = auc_RF>0.5 & p_RF<0.05;
     sel_bad_prefer = auc_RF<0.5& p_RF<0.05;

   %RF

   psth_good = PSTH_TgtOn_Good_RF_stable_Mo{i};
   psth_bad = PSTH_TgtOn_Bad_RF_stable_Mo{i};


   psth_total = [psth_good,psth_bad];
   [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(psth_total');

   PC3 = COEFF(:,1:3);
   psth_good_trans = SCORE(1:size(psth_good,1),1:3);
   psth_bad_trans = SCORE(size(psth_good,1)+1:end,1:3);

   figure
%
   plot3(psth_good_trans(:,1),psth_good_trans(:,2),psth_good_trans(:,3),'-r');
   hold on
   plot3(psth_bad_trans(:,1),psth_bad_trans(:,2),psth_bad_trans(:,3),'-b');
%}
   

%}

  % keyboard

   %{

    psth_goodprefer_good_mean = mean(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_good_prefer,:),1);
    psth_goodprefer_bad_mean = mean(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_good_prefer,:),1);

    psth_goodprefer_good_sem = std(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_good_prefer,:),1)/sqrt(sum(sel));
    psth_goodprefer_bad_sem = std(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_good_prefer,:),1)/sqrt(sum(sel));


     psth_badprefer_good_mean = mean(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_bad_prefer,:),1);
    psth_badprefer_bad_mean = mean(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_bad_prefer,:),1);

    psth_badprefer_good_sem = std(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_bad_prefer,:),1)/sqrt(sum(sel));
    psth_badprefer_bad_sem = std(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_bad_prefer,:),1)/sqrt(sum(sel));


%Non-RF

     psth_goodprefer_good_nonRF_mean = nanmean(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_good_prefer,:),1);
    psth_goodprefer_bad_nonRF_mean = nanmean(PSTH_TgtOn_Bad_nonRF_stable_Mo{i}(sel_good_prefer,:),1);

    psth_goodprefer_good_nonRF_sem = nanstd(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_good_prefer,:),[],1)/sqrt(sum(sel));
    psth_goodprefer_bad_nonRF_sem = nanstd(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_good_prefer,:),[],1)/sqrt(sum(sel));


     psth_badprefer_good_nonRF_mean = nanmean(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_bad_prefer,:),1);
    psth_badprefer_bad_nonRF_mean = nanmean(PSTH_TgtOn_Bad_nonRF_stable_Mo{i}(sel_bad_prefer,:),1);

    psth_badprefer_good_nonRF_sem = nanstd(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_bad_prefer,:),[],1)/sqrt(sum(sel));
    psth_badprefer_bad_nonRF_sem = nanstd(PSTH_TgtOn_Bad_nonRF_stable_Mo{i}(sel_bad_prefer,:),[],1)/sqrt(sum(sel));





     subplot(4,4,(i-1)*4+1);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_good_mean,psth_goodprefer_good_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_bad_mean,psth_goodprefer_bad_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

         subplot(4,4,(i-1)*4+2);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_good_mean,psth_goodprefer_good_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_bad_mean,psth_goodprefer_bad_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);


         subplot(4,4,(i-1)*4+3);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_good_nonRF_mean,psth_goodprefer_good_nonRF_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_bad_nonRF_mean,psth_goodprefer_bad_nonRF_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);



         subplot(4,4,(i-1)*4+4);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_good_nonRF_mean,psth_badprefer_good_nonRF_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_bad_nonRF_mean,psth_badprefer_bad_nonRF_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

   %}



     %{

     good_p = auc_RF(auc_RF>0.5 & p_RF<0.05);
     good_p_non = auc_nonRF(auc_RF>0.5 & p_RF<0.05);

     bad_p = auc_RF(auc_RF<0.5 & p_RF<0.05 );
     bad_p_non = auc_nonRF(auc_RF<0.5 & p_RF<0.05);

     
    plot(auc_nonRF,auc_RF,'ok','LineWidth',1);
    hold on
     plot(good_p_non,good_p,'or','LineWidth',3);
     hold on
     plot(bad_p_non,bad_p,'ob','LineWidth',3);
     

     





 end 
  %}
%{
keyboard;

Group = 2;

good = PSTH_TgtOn_Good_RF_stable_Mo{Group};
bad = PSTH_TgtOn_Bad_RF_stable_Mo{Group};

good_nonRF = PSTH_TgtOn_Good_nonRF_stable_Mo{Group};
bad_nonRF = PSTH_TgtOn_Bad_nonRF_stable_Mo{Group};

goodsac = PSTH_SacOn_Good_RF_stable_Mo{Group};
badsac = PSTH_SacOn_Bad_RF_stable_Mo{Group};

goodsac_nonRF = PSTH_TgtOn_Good_nonRF_stable_Mo{Group};
badsac_nonRF = PSTH_TgtOn_Bad_nonRF_stable_Mo{Group};

 p_Tgt = p_Tgt_RF_stable_Mo{Group};
p_Sac = p_Sac_RF_stable_Mo{Group} ;


name = size(good,1);


for ii = 0:name-1

        fig_index = floor(ii/16)+1;
        subplot_index = mod(ii,16)+1;
        figure(fig_index);
        subplot(4,4,subplot_index);
%
        plot(PSTH_Time_TgtOn_stable,good(ii+1,:),'-r');
        hold on
        plot(PSTH_Time_TgtOn_stable,bad(ii+1,:),'-b');
        plot(PSTH_Time_TgtOn_stable,good_nonRF(ii+1,:),'-m' );
        plot(PSTH_Time_TgtOn_stable,bad_nonRF(ii+1,:),'-g' );

      %{
        plot(PSTH_Time_SacOn_stable,goodsac(ii+1,:),'-r');
        hold on
        plot(PSTH_Time_SacOn_stable,badsac(ii+1,:),'-b');
        plot(PSTH_Time_SacOn_stable,goodsac_nonRF(ii+1,1:120),'-m' );
        plot(PSTH_Time_SacOn_stable,badsac_nonRF(ii+1,1:120),'-g' );
        %}
        if p_Tgt(ii+1)<0.05
            Vis = 1;
        else
            Vis = 0;
        end

        if p_Sac(ii+1)<0.05
            Sac = 1;
        else
            Sac = 0;
        end
        
        


       
       % index=index+1;
      title(sprintf("p_Tgt: %d,p_Sac: %d",Vis,Sac));

 end

%}

%{
PSTH_TgtOn_Good_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_RF_stable(sel)');
            PSTH_TgtOn_Bad_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_RF_stable(sel)');

            PSTH_TgtOn_Good_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_nonRF_stable(sel)');
            PSTH_TgtOn_Bad_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_nonRF_stable(sel)');

%}














%% Plot groups of opto-genetic modulation

if ShowFigureFlag

    %Color order for different groups of neurons

    ColorOrder = [1,0,0;
                  0.9098    0.7020    0.1373;
                  0,0,1;
                  0.3569    0.3451    0.3451];



    FigureStartNum=100;
    FigureIndex=1;

    if PlotOptoEffect

    figtitlestr{FigureIndex}='PSTH_OptoModulation_Type';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    

    for i = 1:length(UniqueMol)
        %subplot(1,length(UniqueMol),i);
        subplot(2,2,i);
        AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

        shadedErrorBar(PSTH_Opto_Time,PSTH_Stim_Mean_Group(i,:),PSTH_Stim_Sem_Group(i,:),'lineprops',{ColorOrder(i,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
        
        legend(sprintf('N=%d', ModuTypeCount(i)));
        
        title(ModulString(i));
        

        set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;
        

    end


    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

    figtitlestr{FigureIndex}='Latency_Dist';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});

    %Latency
    subplot(1,2,1);

    b = bar(Latency_Mean_Group(1:3),'FaceColor','flat');

    b.CData = ColorOrder(1:3,:);

    xticks(UniqueMol(1:3));
    xticklabels({'Fast Ex','Slow Ex','In'});
    
    ylabel('Latency(ms)');


    hold on
    errorbar(UniqueMol(1:3),Latency_Mean_Group(1:3),Latency_Sem_Group(1:3),'ok','LineWidth',3);
    
    box off
    set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

    %Strength
    subplot(1,2,2);

    b = bar(MS_Mean_Group(1:3),'FaceColor','flat');

    b.CData = ColorOrder(1:3,:);

    xticks(UniqueMol(1:3));
    xticklabels({'Fast Ex','Slow Ex','In'});
    
    ylabel('Strength');


    hold on
    errorbar(UniqueMol(1:3),MS_Mean_Group(1:3),MS_Sem_Group(1:3),'ok','LineWidth',3);
    
    box off
    set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');


%% NonSigGroup
 FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

  figtitlestr{FigureIndex}='PSTH_NonSigSep';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,400],'Name',figtitlestr{FigureIndex});

   ColorOrder_NonSig=[1,0,0;0,0,1;0.5,0,0;0,0,0.5];

    %PSTH_NonSigStim_Mean_Group

    subplot(1,2,1);

     AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(1,:),PSTH_NonSigStim_Sem_Group(1,:),'lineprops',{ColorOrder_NonSig(1,:)},'transparent',1,'patchSaturation',0.3);
        hold on
        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(2,:),PSTH_NonSigStim_Sem_Group(2,:),'lineprops',{ColorOrder_NonSig(2,:)},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        
        %legend(sprintf('N=%d', NonSigCount(i)));
        legend({sprintf('Ex;N=%d',NonSigCount(1)),sprintf('In;N=%d',NonSigCount(2))});
        
       % title(ModulString(i));
        

        set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;


    subplot(1,2,2);

       AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(3,:),PSTH_NonSigStim_Sem_Group(3,:),'lineprops',{ColorOrder_NonSig(3,:)},'transparent',1,'patchSaturation',0.3);
        hold on
        try
        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(4,:),PSTH_NonSigStim_Sem_Group(4,:),'lineprops',{ColorOrder_NonSig(4,:)},'transparent',1,'patchSaturation',0.3);
        catch
            disp('No enough data to show group 4')
        end
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        
        %legend(sprintf('N=%d', NonSigCount(i)));
        try
        legend({sprintf('Ex;N=%d',NonSigCount(3)),sprintf('In;N=%d',NonSigCount(4))});
        catch
            legend({sprintf('Ex;N=%d',NonSigCount(3))});
        end
        
       % title(ModulString(i));
        

        set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;

        %% Latency Distribution
         FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

        figtitlestr{FigureIndex}='Latency_Distribution';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
   Edges =[0:5:200];

   subplot(2,1,1);
   histogram(Latency_Ex,Edges,'Facecolor','r','LineWidth',3);
   hold on
   histogram(Latency_Ex(2:end),Edges(2:end),'Facecolor','y','LineWidth',3);

  % set(gca, 'YScale', 'log');
   xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   box off



   subplot(2,1,2);
   histogram(Latency_In,Edges,'Facecolor','b','LineWidth',3);
   %set(gca, 'YScale', 'log');
   xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    box off


%% Pie plot 
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

        figtitlestr{FigureIndex}='Modulation';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});

   h=piechart([FastExcitationNum,SlowExcitationNum,InhibitionNum,NonSigNum],["Fast Excitation","Slow Excitation","Inhibition","No-Effect"]);
   h.ColorOrder = ColorOrder;
%{

    ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);

%}
%% Modulation Strength
%{
MS_Mean_Group(i) = mean(abs(SM(sel)));
    MS_Sem_Group(i) = std(abs(SM(sel)))/sqrt(sum(sel));
%}

FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;

 figtitlestr{FigureIndex}='ModulationStrength';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
        
   b = bar([1,2,3],MS_Mean_Group(1:3));



    
   
    end %End of if plot opto effect



end %End of if figureplot
end %End of the function

  function checkGroup(name,data)
        for ii = 0:length(name )-1

        fig_index = floor(ii/16)+1;
        subplot_index = mod(ii,16)+1;
        figure(fig_index);
        subplot(4,4,subplot_index);
    

       plot(data(ii+1,:));
       % index=index+1;
       title(name(ii+1))

        end


    end

 