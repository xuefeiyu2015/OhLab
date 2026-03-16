
function Batch_OptogeneticsFreeView_SC(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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

 Monkey = FilesName_Each{i}(1);

       if Monkey == 'R'
           MonkeyIndex(index) = 1;%Monkey R
           MonkeyIndex_session(i) = 1;
       else
           MonkeyIndex(index) = 2;%Monkey A
           MonkeyIndex_session(i) = 2;
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
%     p_AS(index) = Data('Ac_p');
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



     %Memory saccade
     if isfield(OutputData,'MemorySaccade')
     
     MemorySaccadeIndex(index) = 1;
     
     Data =OutputData.MemorySaccade(1).DataStamp;

     %PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF');
     %PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF');

     Mean_FR_MS(index) = mean([nanmean(Data('PSTH_TgtOn_RF'),'all'),nanmean(Data('PSTH_SacOn_RF'),'all')]);
     

     PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF_z');
     PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF_z');

     try
    RF_h_v(index) = Data('RF_h_v');
    RF_h_m(index) = Data('RF_h_m');
catch
    tmp = Data('RF_h_v');
    RF_h_v(index) =tmp(1);
    tmp = Data('RF_h_m');
    RF_h_m(index)=tmp(1);
     end


     %CellType(index) = Data('CellType');%1:visual;2:motor;3:vis-motor
     %Reasign cell_type
     %
     if RF_h_v(index)==1 & (RF_h_m(index)~=1)

         CellType(index) = 1;%Visual-Only

     elseif RF_h_v(index)==1 & RF_h_m(index)==1
         CellType(index) = 2;% Visual-Motor

     elseif (RF_h_v(index)~=1) & RF_h_m(index)==1
         CellType(index) = 3;% Motor-Only

     elseif RF_h_v(index)==0 & RF_h_m(index)==0

         CellType(index) = 5;%Null Type 
     else
         CellType(index) = 4; %Other type
     end


     


     

PSTH_Time_TgtOn = Data('PSTH_Time_TgtOn');
     PSTH_Time_SacOn = Data('PSTH_Time_SacOn');
     else

         MemorySaccadeIndex(index) = 0;

     PSTH_Tgt_RF{index} = NaN*ones(1,140);
     PSTH_Sac_RF{index} = NaN*ones(1,140);


     

     RF_h_v(index) = NaN;
     RF_h_m(index) = NaN;

     CellType(index) =NaN;%1:visual;2:motor;3:vis-motor

   

     
         Mean_FR_MS(index)=NaN;
     end


 
     index = index+1;   
    end %End of for file name
     
end


%Plot controller
PlotOptoEffect = 1;
PlotMemorySaccade =1;


%Analysis controller 
AnalysisOS = 0;% For Video Free View
AnalysisMS = 1;% For memory saccade


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

MonkeyIndex=MonkeyIndex(ScreenCell);


%% For VideoFreeView

N = length(MonkeyIndex);

N_R = sum(MonkeyIndex == 1);
N_A = sum(MonkeyIndex == 2);

disp('Total Number: ')
disp(N)
disp('Monkey R total number:')
disp(N_R)
disp('Monkey A total number:')
disp(N_A)

disp('Sessions')
disp('Monkey R total session: ')
disp(sum(MonkeyIndex_session==1))
disp('Monkey A total session: ')
disp(sum(MonkeyIndex_session==2))





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
p_Cri50=p_Cri50(ScreenCell);

%maxnum_ms = max(cellfun(@numel,PSTH_Tgt_RF));


PSTH_Tgt_RF = cell2mat(cellfun(@(x) x(1:140), PSTH_Tgt_RF(ScreenCell),'UniformOutput',false)');
PSTH_Sac_RF = cell2mat(PSTH_Sac_RF(ScreenCell)');

PSTH_Time_TgtOn = PSTH_Time_TgtOn(1:140);

MemorySaccadeIndex = MemorySaccadeIndex(ScreenCell)';
%p_TgtV_prefer=p_TgtV_prefer(ScreenCell,:);

%p_SacV_prefer =  p_SacV_prefer(ScreenCell,:);



%VisTuning_Mean = VisTuning_Mean(ScreenCell,:);

%VisTuning_Sem = VisTuning_Sem(ScreenCell,:);



%SacTuning_Mean = SacTuning_Mean(ScreenCell,:);

%SacTuning_Sem = SacTuning_Sem(ScreenCell,:);

Mean_FR_MS=Mean_FR_MS(ScreenCell)';
CellType = CellType(ScreenCell)' ;







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


Fast_M = Latency<10;
Slow_M = Latency>=10;
Excitation = SM> 0;
Inhibition = SM<=0;

%Sig_m = SM_p < 0.05;
%Sig_m = ~isnan(Latency) & AS >=5;

%Sig_m = (p_whole_cri50<0.05 | SM_p<0.05 | PMS_p<0.05);%&(~isnan(Latency));
Sig_m = SM_p<0.05 & (~isnan(Latency));
%Sig_m = p_whole_cri50>-0.1;
%Latency
Non_Sig_m = ~Sig_m;

SepLevel = 0.1;

Non_Sig_level2 = p_Cri50>SepLevel & SM_p>SepLevel;
Non_Sig_level1 = ~Non_Sig_level2  & Non_Sig_m;

%Non_Sig_m = Non_Sig_level2;


Ex_Fast_Sig = Fast_M & Excitation & Sig_m;

Ex_Slow_Sig = Slow_M & Excitation & Sig_m;

Ex_Sig = Ex_Fast_Sig | Ex_Slow_Sig;


%In_Fast_Sig = Fast_M & Inhibition & Sig_m;

%In_Slow_Sig = Slow_M & Inhibition & Sig_m;

In_Sig = Inhibition & Sig_m;

PSTH_Use = PSTH_norm;
%PSTH_Use = PSTH; 
%PSTH_Use = PSTH_z; 
%PSTH_Use = PSTH_norm2;

ModuType = NaN*ones(1,length(Sig_m));
NSig_Level = NaN*ones(1,length(Sig_m));



%%
SeparateGroup=4;

if SeparateGroup == 3 %Ex,Inh,non

    ModuType(Ex_Sig) =1;
    %ModuType(Ex_Slow_Sig) =2;

    ModuType(In_Sig) =3;
    %ModuType(In_Slow_Sig) =4;

    ModuType(Non_Sig_m) =4;
    ModulString = {'Excitation','Inhibition','No Effect'};

     ColorOrder = [1,0,0;
                %  0.9098    0.7020    0.1373;
                  0,0,1;
                  0.3569    0.3451    0.3451];

else  %Separate into 4

    ModuType(Ex_Fast_Sig) =1;
    ModuType(Ex_Slow_Sig) =2;

ModuType(In_Sig) =3;
%ModuType(In_Slow_Sig) =4;

ModuType(Non_Sig_m) =4;
ModulString = {'Fast Excitation','Slow Excitation','Inhibition','No Effect'};
 ColorOrder = [1,0,0;
                  0.9098    0.7020    0.1373;
                  0,0,1;
                  0.3569    0.3451    0.3451];

end

NSig_Level(Non_Sig_level1&Excitation) =1;
NSig_Level(Non_Sig_level1&Inhibition) =2;

NSig_Level(Non_Sig_level2&Excitation) =3;
NSig_Level(Non_Sig_level2&Inhibition) =4;

NSig_Level=NSig_Level';


%ModulString = {'Fast Excitation','Slow Excitation','Inhibition','No Effect'};

%NonSigString = {'NonSig_Ex_l1','NonSig_In_l1','NonSig_Ex_l2','NonSig_In_l2'};


NonSigString = {'NonSig_Ex_l1','NonSig_In_l1','NonSig_Ex_l2','NonSig_In_l2'};

Latency_Ex = Latency(Ex_Fast_Sig|Ex_Slow_Sig);
Latency_In = Latency(In_Sig);



UniqueMol = unique(ModuType);

UniqueMol = UniqueMol(~isnan(UniqueMol));

ModuType=ModuType';%Transfer to column;

UniqueNonSig = unique(NSig_Level);
UniqueNonSig=UniqueNonSig(~isnan(UniqueNonSig));

%For bar plot
ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);

%

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
%cla reset; % Clear the current axes and reset properties
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




if AnalysisMS

    MemorySaccade_Sel = MemorySaccadeIndex ==1;

    CellTypeUnique = unique(CellType);
    CellTypeUnique = CellTypeUnique(~isnan(CellTypeUnique));

    %CellTypeUnique = CellTypeUnique(~isnan(CellTypeUnique) & CellTypeUnique<=3);

   % CellTypeUnique = ones(size(CellTypeUnique));%For a test

  % A=[Mean_FR_MS',MemorySaccadeIndex,CellType];

  Select_FR = Mean_FR_MS > 1;
  disp(sprintf('Memory Saccade exclude: %1.2f',sum(Select_FR==0)/length(Select_FR)));

  % UniqueMol = 1;%For a test;
 % ModuType = ones(size(ModuType));%For a test
    for i = 1:length(UniqueMol)
        for j = 1:length(CellTypeUnique)
           %x
           % if UniqueMol(i) == 1 &  (CellTypeUnique(j)==3)
           % sel = (ModuType == 1| ModuType == 2) & CellType == CellTypeUnique(j) & MemorySaccade_Sel & Select_FR;
           % sel = (ModuType == 1| ModuType == 2) & CellType == CellTypeUnique(j) & MemorySaccade_Sel & Select_FR;
          %sel = ModuType == UniqueMol(i) & CellType == CellTypeUnique(j) & MemorySaccade_Sel & Select_FR;
          
          %  else
            sel = ModuType == UniqueMol(i) & CellType == CellTypeUnique(j) & MemorySaccade_Sel & Select_FR;
         %   end
           %}
           % sel = ModuType == UniqueMol(i) & CellType == CellTypeUnique(j) & MemorySaccade_Sel & Select_FR;

            TypeCounts(i,j) = sum(sel);

            PSTH_Tgt_Mean{i,j} = nanmean(PSTH_Tgt_RF(sel,:),1);
            PSTH_Sac_Mean{i,j} = nanmean(PSTH_Sac_RF(sel,:),1);

            PSTH_Tgt_Sem{i,j} = nanstd(PSTH_Tgt_RF(sel,:),[],1)/sqrt(sum(sel));
            PSTH_Sac_Sem{i,j} = nanstd(PSTH_Sac_RF(sel,:),[],1)/sqrt(sum(sel));

           % Mean_FR_MS


%Copy from here
           PSTH_Tgt_Each{i,j} = PSTH_Tgt_RF(sel,:);
            PSTH_Sac_Each{i,j} = PSTH_Sac_RF(sel,:);

           
        

        end
     end


%devider = repmat(sum( TypeCounts,2),1,size(TypeCounts,2));

devider = repmat(sum( TypeCounts,1),size(TypeCounts,1),1);
    
 TypeFreq =  TypeCounts./devider;

 %Statistics:
 Class1_total = devider(1,1:3);
 Class1_group =TypeCounts(1,1:3);
 Class1_others =Class1_total -Class1_group;
 observed = [Class1_group; Class1_others];

 

% Use built-in function if available
[p_ms_class1,Q] = chi2test(observed);



keyboard
 
end 




%}
%{
%% For nonsig neurons 
for i = 1:length(UniqueNonSig)
    sel = NSig_Level == UniqueNonSig(i);
    NonSigCount(i) = sum(sel);


    PSTH_NonSigStim_Group = PSTH_Use(sel,:);

   

  

    PSTH_NonSigStim_Mean_Group(i,:) = mean(PSTH_NonSigStim_Group,1);

    PSTH_NonSigStim_Sem_Group(i,:) = std(PSTH_NonSigStim_Group,[],1)/sqrt(sum(sel));

   



end

%}












%% Plot groups of opto-genetic modulation

if ShowFigureFlag

    %Color order for different groups of neurons

   



    FigureStartNum=100;
    FigureIndex=1;

    if PlotOptoEffect
if SeparateGroup == 3
    figtitlestr{FigureIndex}='PSTH_OptoModulation_Excitation';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [2 2], 'nature');

    subplot(2,1,1);
       AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        %{
        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        %}
        shadedErrorBar(PSTH_Opto_Time,PSTH_Stim_Mean_Group(UniqueMol==1,:),PSTH_Stim_Sem_Group(UniqueMol==1,:),'lineprops',{ColorOrder(UniqueMol==1,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);
        ylim([-0.1,1]);
       % legend(sprintf('N=%d', ModuTypeCount(1)));
        
       % title('Population');
       % xlabel('Time from Opto-stim onset (ms)')
        ylabel('Normalized Firing Rate(Hz)');
        

        set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
      %  xlabel('Time from Opto-stim Onset (ms)');
      %  ylabel('Normalized Responses');

        box off;

%{
%% Inhibition population
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;
         figtitlestr{FigureIndex}='PSTH_OptoModulation_Inhibition';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
%}
        subplot(2,1,2);
       AverageRegion = [0,100];
        area(AverageRegion,[-0.3,-0.3],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        %
        AverageRegion = [0,100];
        area(AverageRegion,[0.2,0.2],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        %
        shadedErrorBar(PSTH_Opto_Time,PSTH_Stim_Mean_Group(UniqueMol==3,:),PSTH_Stim_Sem_Group(UniqueMol==3,:),'lineprops',{ColorOrder(UniqueMol==3,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);
        ylim([-0.3,0.3]);
       % legend(sprintf('N=%d', ModuTypeCount(1)));
        
       % title('Population');
        xlabel('Time from Opto-stim onset (ms)')
        ylabel('Normalized Firing Rate(Hz)');
        

        set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
        xlabel('Time from Opto-stim Onset (ms)');
        ylabel('Normalized Responses');

        box off;
else



    figtitlestr{FigureIndex}='PSTH_OptoModulation_Type';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 2], 'nature');
    AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');


    for i = 1:length(UniqueMol)
        %subplot(1,length(UniqueMol),i);
       % subplot(2,2,i);
      

        shadedErrorBar(PSTH_Opto_Time,PSTH_Stim_Mean_Group(i,:),PSTH_Stim_Sem_Group(i,:),'lineprops',{ColorOrder(i,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);
        
       
        
      %  title(ModulString(i));
        

      
        % legstring{i}=sprintf('N=%d', ModuTypeCount(i));

    end

   

     legend({'Fast Ex','Slow Ex','Inh','No effect'});
       set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;


end

        
%keyboard
%{
    for i = 1:length(UniqueMol(UniqueMol<4))
        %subplot(1,length(UniqueMol),i);
        subplot(2,1,i);
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
%}
%{
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
%}
%{
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
        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(4,:),PSTH_NonSigStim_Sem_Group(4,:),'lineprops',{ColorOrder_NonSig(4,:)},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        
        %legend(sprintf('N=%d', NonSigCount(i)));
        legend({sprintf('Ex;N=%d',NonSigCount(3)),sprintf('In;N=%d',NonSigCount(4))});
        
       % title(ModulString(i));
        

        set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;
%}
    %{    
    %% Latency Distribution
         FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

        figtitlestr{FigureIndex}='Latency_Distribution';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
   Edges =[0:5:200];

   subplot(2,1,1);
   histogram(Latency_Ex,Edges,'Facecolor','r','LineWidth',3);
   set(gca, 'YScale', 'log');
   xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   box off



   subplot(2,1,2);
   histogram(Latency_In,Edges,'Facecolor','b','LineWidth',3);
   set(gca, 'YScale', 'log');
   xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    box off
    %}

%% Bar plot 
if SeparateGroup == 3
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

        figtitlestr{FigureIndex}='Modulation';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [2 2], 'nature');

   data = [ExcitationNum,InhibitionNum,NonSigNum];
   %h=piechart(data,["Excitation","Inhibition","No-Effect"]);
   h = bar([1,2,3],data);
    h.FaceColor="flat";
 h.CData=[1,0,0;0,0,1;0,0,0];
 set(gca,'XTickLabel',["Excitation","Inhibition","No-Effect"]);
 ylabel('Number of neurons');
 set(gca,'LineWidth',0.5,'FontSize',5);

 box off;
else
    %% Pie plot 
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

   


        figtitlestr{FigureIndex}='Modulation';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 300,300],'Name',figtitlestr{FigureIndex});
% set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
  %  formatFig(fig{FigureIndex}, [0.5 0.5], 'nature');

   h = piechart(ModuTypeCount,["FastExc","SlowExc","Inh","No-Effect"])%,'FontSize',5);%,["FastExc","SlowExc","Inh","No-Effect"]);
   %
   %{
  for i = 1:length(ModuTypeCount)
      modu_label{i} = sprintf('N = %d',ModuTypeCount(i));
  end
  h.Labels = modu_label;
   %}
   %h.Labels = '';
  h.ColorOrder = ColorOrder;
   h.StartAngle = -30;
end

   %keyboard


   % Apply the colors to each segment
 %  h.ColorOrder = ColorOrder;


%{

    ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);

%}
%% Modulation Strength
%{
MS_Mean_Group(i) = mean(abs(SM(sel)));
    MS_Sem_Group(i) = std(abs(SM(sel)))/sqrt(sum(sel));


FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;

 figtitlestr{FigureIndex}='ModulationStrength';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
     
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');

   b = bar([1,2,3],MS_Mean_Group(1:3));
%}



    
   
    end %End of if plot opto effect


    %% Memory saccade
  if  PlotMemorySaccade 
%% Population psth

 %for i = 1: length(UniqueMol)
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

    figtitlestr{FigureIndex}= sprintf('PSTH_MemorySac_Type:%d',UniqueMol(1));
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
     set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 3.5], 'nature');

    typestring = {'Visual','Visual-Movement','Movement','Others'};

   for j = 1:length(CellTypeUnique)-1
        subplot(length(CellTypeUnique)-1,2,2*(j-1)+1);
      % subplot(1,2,1);

        %% Target onset 


        shadedErrorBar(PSTH_Time_TgtOn ,PSTH_Tgt_Mean{1,j},PSTH_Tgt_Sem{1,j},'lineprops',{[128,128,128]/255},'transparent',1,'patchSaturation',0.2)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);
   %      ylim([-0.5,3]);
       xlim([-250,300]);
       hold on
       
       ymax = ceil(max(max(PSTH_Tgt_Mean{1,j}+PSTH_Tgt_Sem{1,j}),max(PSTH_Sac_Mean{1,j}+PSTH_Sac_Sem{1,j})));
       ymin = floor(min(min(PSTH_Tgt_Mean{1,j}-PSTH_Tgt_Sem{1,j}),min(PSTH_Sac_Mean{1,j}-PSTH_Sac_Sem{1,j})));
       ylim([ymin,ymax]);
       plot([0,0],[ymin,ymax],'--k','LineWidth',0.5);
       if j == 2
            ylabel('Z-Score');
        end
        if j == 3
            xlabel('Time from event onset (ms)');
        end
        title(typestring(j));
        set(gca,'FontSize',5);
        

        %% Saccade onset
        subplot(length(CellTypeUnique)-1,2,2*(j-1)+2);
 % subplot(1,2,2);

        %Saccade on 

        shadedErrorBar(PSTH_Time_SacOn ,PSTH_Sac_Mean{1,j},PSTH_Sac_Sem{1,j},'lineprops',{[128,128,128]/255},'transparent',1,'patchSaturation',0.2)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);

         xlim([-450,200]);
        ylim([ymin,ymax]);
        hold on
        plot([0,0],[ymin,ymax],'--k','LineWidth',0.5);
        set(gca,'FontSize',5);
      %  ylim([-0.5,3]);


%PSTH_Tgt_Mean{i,j}


   end
   %Bar plot for fast activated neuron only
 FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

     figtitlestr{FigureIndex}='PSTH_MemorySac_ProjectionNeuron';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 3.5], 'nature');
    subplot(2,1,1);
     b = bar(TypeCounts(1,1:3),'FaceColor','k');

     
%         b.FaceColor = 'Flat';
        % b.CData = MS_TypeOrder(1,:);
     

     box off
     xticks(UniqueMol);
    xticklabels({'Vis','Vis-Mov','Mov'});
    ylabel('Number of neurons');
    
  %  legend({'Vis-Only','Vis-Motor','Motor'});

     set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');

subplot(2,1,2);
%% Frequency
      b = bar(TypeFreq(1,1:3),'FaceColor','k');

     
%         b.FaceColor = 'Flat';
        % b.CData = MS_TypeOrder(1,:);
     

     box off
     xticks(UniqueMol);
    xticklabels({'Vis','Vis-Mov','Mov'});
    ylabel('Proportion');
  %  legend({'Vis-Only','Vis-Motor','Motor'});

     set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
    % ylim([0,0.6]);
 end %End of the if plot memory saccade



end %End of if figureplot

%% Output to files
%Setup output directory
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
%MarkerFolder='DataHub';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));


%%%%%%%%%%%Export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up output path
%OutputFolerName='ForagingTask';
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end),'ExportFigure',BasicDirectory(end));%,OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
end
cd(OutPath);
%{
FileName=Data.FileName;
FileName=FileName(1:end-4);


%Get File Information
num = regexp(FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));

NeuronNum=cell2mat(num(2));
ChannelNum=SpikeChannel;


%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);
%}

%Output file name
%OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
%OutputFigureName=sprintf('%s_%s_N%s_C%s_%CID_%NID',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),string(ChannelID),string(ClusterID));

OutputFigureName='SC_VideoVeiw_Neural_V2';

for jj=1:FigureIndex
        if ~isempty(fig{jj})
        

         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.pdf');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
end

    disp('Figures have been exported to the exported folder');
%keyboard
%{
ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'FreeViewNeural')
    ExistFlag=1;
    end
end



%{
%FileName=strcat(FileName,'_RearchingTime.mat');
if SaveFlag
save(FileName,'OutputData');
disp('Seaching Time Data Exported');
else
    disp('No data export');
end
%}
%OutputData For current analysis
OutputData_New=[];

OutputData_New.FreeViewNeural.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data


NameStr={'PSTH_Time','PSTH_Stim_Mean','PSTH_Control_Mean','AcStrength','ResponseLatency',...
    'PostModulationStrength','PostModulation_p','StimModulation','StimModulation_p',...
    'PSTH_FirstHalf','PSTH_SecondHalf',...
    'AS_First','AS_Second',...
    'PSTH_Stim_Mean_z','PSTH_Control_Mean_z','PSTH_Stim_Mean_First_z','PSTH_Stim_Mean_Second_z',...
    'StimModulation_z',...
    'PSTH_Stim_Mean_norm','PSTH_Control_Mean_norm','FR_Mean','p_StimControl',...
    'AcDur','psth_stim_normed2','psth_control_normed2','p_cri50','p_cri20','Strength_sig','p_whole_cri50',...
    'Latency3','Latency2'...
    'Raster_Stim','Raster_Control'
    };

DataLib={MeanPSTH_Time_Mean,SpikeCountStim_Mean,SpikeCountControl_Mean,ActivationStrength,ResponseLatency,...
    PostStimModulation,p_PostModulation,StimModulation,p_StimModulation,...
    PSTH_Stim_Mean_First,PSTH_Stim_Mean_Second,...
    ActivationStrength_FirstHalf,ActivationStrength_SecondHalf,...
    PSTH_Stim_Mean_z,PSTH_Control_Mean_z,PSTH_Stim_Mean_First_z,PSTH_Stim_Mean_Second_z,...
    StimModulation_z,...
    PSTH_Stim_Mean_norm_raw,PSTH_Control_Mean_norm_raw,FR_Mean,p_StimControl,...
    AcDur,psth_stim_normed2,psth_control_normed2,p_use_cri,p_use_cri20,Strength_sig,p_whole_cri50,...
    Latency3,Latency2,...
    Raster_Stim,Raster_Control



};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.FreeViewNeural.Task=TaskType;
OutputData_New.FreeViewNeural.TaskCode=TaskCode;
OutputData_New.FreeViewNeural.DataStamp=DataStamp;

OutputData_New.FreeViewNeural.StartTrial =StartTrial;
%{
%Store figures
if ShowFigureFlag
   
     %OutputData_New.SceneTuning.Figures={fig99,fig100,fig101,fig102};
    %Output figures in a folder;
    FolderName=OutputFigureName;
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
    
    TotalFigure=length(fig);
    
    for jj=1:TotalFigure
        if ~isempty(fig{jj})
        
         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.jpg');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
        
    end

    disp('Figures have been exported to the neuron folder');

end

%}





if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.FreeViewNeural.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
     elseif ismember(FileName,{'Adams062124VIDEOVIEW1'})
         %May save different task segment
         if StartTrial~=OutputData.FreeViewNeural.StartTrial

             if StartTrial > 160
                 OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum+1),string(ChannelNum));
                 OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             elseif StartTrial > 360
                  OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum+2),string(ChannelNum));
                  OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             end


         end

          
     else
         %If the same task, replace the old one with the new one
         %OutputData=OutputData_New;
         OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;
         
    
    
     end
else
    %Output the current file
   OutputData.FreeViewNeural=OutputData_New.FreeViewNeural; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  
%}
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

 