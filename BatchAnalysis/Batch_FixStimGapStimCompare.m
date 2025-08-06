function Batch_FixStimGapStimCompare(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
FilesName_Each=Data.FileName;

UseAll = 0;%1: Use all the recorded sites; 0: only use significant sites 

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
  
     load(FilesName{i}{1});
     
     OutputDataTemp=OutputData;
     
     if isfield(OutputDataTemp,'FreeViewBehavior') && isfield(OutputDataTemp,'FixViewBehavior') &&  isfield(OutputDataTemp,'GapFixViewBehavior')
         
     
    
    
%

      Monkey = FilesName_Each{i}(1);

          if Monkey == 'R'
              MonkeyIndex(index) = 1;
          else
              MonkeyIndex(index) = 2;
          end

     %% Free View
     %Task(i)=OutputData.FreeViewBehavior(1).TaskCode;
     %TaskName{i}=OutputData.FreeViewBehavior(1).Task;
     
     
     Data =OutputData.FreeViewBehavior.DataStamp;
     p_FreeView(index) = Data('p_prop');

    % PropClSac(i) = Data('PropClSac');
    % p_prop(i) = Data('p_prop');

     NumberOfAngle_Control(index,:) = Data('NumberOfAngle_Control');
     NumberOfAngle_Stim(index,:) = Data('NumberOfAngle_Stim');


     tmp = Data('AngleAmpEachTrial');

     Angle_each{index} = tmp(:,1)';
     Amp_each{index} = tmp(:,2)';

     tmp = Data('AngleAmpEachTrial_Control');

     Angle_each_Control{index} = tmp(:,1)';
     Amp_each_Control{index} = tmp(:,2)';

     
%{
     SaccadeCountStim_Sum(index,:) = Data('SaccadeCountStim');
     SaccadeCountControl_Sum(index,:) = Data('SaccadeCountControl');
%}

%{
     theta_stim{i}=Data('theta_Stim');
     amp_stim{i}=Data('Amp_Stim');
     theta_control{i}=Data('theta_Control');
     amp_control{i}=Data('Amp_Control');

%}




  %{   
     
     TimeSeq = Data('TimeForSaccadeRaster');

     Data =OutputData.FreeViewNeural.DataStamp;


     SM_p(i) = Data('StimModulation_p');
  %}

 %% FixStim
     %Task(i)=OutputData.FreeViewBehavior(1).TaskCode;
     %TaskName{i}=OutputData.FreeViewBehavior(1).Task;
     
     
     Data =OutputData.FixViewBehavior.DataStamp;
     p_FixStim(index) = Data('p_prop');

    % PropClSac(i) = Data('PropClSac');
    % p_prop(i) = Data('p_prop');

     NumberOfAngle_Control_Fix(index,:) = Data('NumberOfAngle_Control');
     NumberOfAngle_Stim_Fix(index,:) = Data('NumberOfAngle_Stim');

     tmp = Data('AngleAmpEachTrial');

     Angle_each_Fix{index} = tmp(:,1)';
     Amp_each_Fix{index} = tmp(:,2)';

     tmp = Data('AngleAmpEachTrial_Control');

     Angle_each_FixControl{index} = tmp(:,1)';
     Amp_each_FixControl{index} = tmp(:,2)';



     %% GapFixStim
     %Task(i)=OutputData.FreeViewBehavior(1).TaskCode;
     %TaskName{i}=OutputData.FreeViewBehavior(1).Task;
     
     
     Data =OutputData.GapFixViewBehavior.DataStamp;
     p_GapFixStim(index) = Data('p_prop');

    % PropClSac(i) = Data('PropClSac');
    % p_prop(i) = Data('p_prop');

     NumberOfAngle_Control_GapFix(index,:) = Data('NumberOfAngle_Control');
     NumberOfAngle_Stim_GapFix(index,:) = Data('NumberOfAngle_Stim');

      tmp = Data('AngleAmpEachTrial');

     Angle_each_GapFix{index} = tmp(:,1)';
     Amp_each_GapFix{index} = tmp(:,2)';


     tmp = Data('AngleAmpEachTrial_Control');

     Angle_each_GapFixControl{index} = tmp(:,1)';
     Amp_each_GapFixControl{index} = tmp(:,2)';

     index=index+1;

     
     end

     end
     %%
 
%end
%{
NumAngleStim = sum(NumberOfAngle_Stim,1);
PropAngleStim = sum(NumberOfAngle_Stim,1)/sum(NumberOfAngle_Stim,'all')*100;

NumAngleStim_Fix = sum(NumberOfAngle_Stim_Fix,1);
PropAngleStim_Fix = sum(NumberOfAngle_Stim_Fix,1)/sum(NumberOfAngle_Stim_Fix,'all')*100;

NumAngleStim_GapFix = sum(NumberOfAngle_Stim_GapFix,1);
PropAngleStim_GapFix = sum(NumberOfAngle_Stim_GapFix,1)/sum(NumberOfAngle_Stim_GapFix,'all')*100;
%}

sig = p_FreeView<0.05 & p_FixStim<0.05 & p_GapFixStim<0.05;


%disp(sum(sig));

Angle_each=Angle_each(sig);
Angle_each_Fix=Angle_each_Fix(sig);
Angle_each_GapFix=Angle_each_GapFix(sig);

Angle_each_Control=Angle_each_Control(sig);
Angle_each_FixControl=Angle_each_FixControl(sig);
Angle_each_GapFixControl=Angle_each_GapFixControl(sig);


Amp_each=Amp_each(sig);
Amp_each_Fix=Amp_each_Fix(sig);
Amp_each_GapFix=Amp_each_GapFix(sig);


Amp_each_Control = Amp_each_Control(sig);
Amp_each_FixControl=Amp_each_FixControl(sig);
Amp_each_GapFixControl=Amp_each_GapFixControl(sig);

 Angle_each_line = cell2mat(Angle_each);
 Angle_each_Fix_line = cell2mat(Angle_each_Fix);
 Angle_each_GapFix_line = cell2mat(Angle_each_GapFix);

  Angle_each_Controlline = cell2mat(Angle_each_Control);
 Angle_each_Fix_Controlline = cell2mat(Angle_each_FixControl);
 Angle_each_GapFix_Controlline = cell2mat(Angle_each_GapFixControl);

 Angle_each_line_valid = Angle_each_line(~isnan(Angle_each_line));
 Angle_each_Fixline_valid = Angle_each_Fix_line(~isnan(Angle_each_Fix_line));
 Angle_each_GapFixline_valid = Angle_each_GapFix_line(~isnan(Angle_each_GapFix_line));

 NumberOfSaccades_FreeView = length(Angle_each_line_valid);
 NumberOfSaccades_Fix = length(Angle_each_Fixline_valid);
 NumberOfSaccades_GapFix = length(Angle_each_GapFixline_valid);
 



Amp_each_line = cell2mat(Amp_each);
 Amp_each_Fix_line = cell2mat(Amp_each_Fix);
 Amp_each_GapFix_line = cell2mat(Amp_each_GapFix);

 Amp_each_line_valid= Amp_each_line(~isnan(Amp_each_line));
 Amp_each_Fixline_valid = Amp_each_Fix_line(~isnan(Amp_each_Fix_line));
 Amp_each_GapFixline_valid = Amp_each_GapFix_line(~isnan(Amp_each_GapFix_line));

 Amp_free_mean = mean(Amp_each_line_valid);
 Amp_free_sem = std(Amp_each_line_valid)/sqrt(length(Amp_each_line_valid));

Amp_fix_mean = mean(Amp_each_Fixline_valid);
 Amp_fix_sem = std(Amp_each_Fixline_valid)/sqrt(length(Amp_each_Fixline_valid));

 Amp_gapfix_mean = mean(Amp_each_GapFixline_valid);
 Amp_gapfix_sem = std(Amp_each_GapFixline_valid)/sqrt(length(Amp_each_GapFixline_valid));

%Statistics for the amplitude

index_task = [ones(1,length(Amp_each_line_valid)),2*ones(1,length(Amp_each_GapFixline_valid)),3*ones(1,length(Amp_each_Fixline_valid))]';
amp_task =[Amp_each_line_valid,Amp_each_Fixline_valid,Amp_each_GapFixline_valid]';

[p_Amp,~,~] = anova1(amp_task,index_task,'off');

 Amp_each_Controlline = cell2mat(Amp_each_Control);
 Amp_each_Fix_Controlline = cell2mat(Amp_each_FixControl);
 Amp_each_GapFix_Controlline = cell2mat(Amp_each_GapFixControl);


 Prop_Sac = (length(Angle_each_line)-sum(isnan(Angle_each_line)))/length(Angle_each_line);
 Prop_Sac_Fix = (length(Angle_each_Fix_line)-sum(isnan(Angle_each_Fix_line)))/length(Angle_each_Fix_line);
 Prop_Sac_GapFix = (length(Angle_each_GapFix_line)-sum(isnan(Angle_each_GapFix_line)))/length(Angle_each_GapFix_line);


 Prop_Sac_Control = (length(Angle_each_Controlline)-sum(isnan(Angle_each_Controlline)))/length(Angle_each_Controlline);
 Prop_Sac_Fix_Control = (length(Angle_each_Fix_Controlline)-sum(isnan(Angle_each_Fix_Controlline)))/length(Angle_each_Fix_Controlline);
 Prop_Sac_GapFix_Control = (length(Angle_each_GapFix_Controlline)-sum(isnan(Angle_each_GapFix_Controlline)))/length(Angle_each_GapFix_Controlline);

%Statistics:
%Free View
Frewview_Stim = [(length(Angle_each_line)-sum(isnan(Angle_each_line))),sum(isnan(Angle_each_line))];
Frewview_Control = [(length(Angle_each_Controlline)-sum(isnan(Angle_each_Controlline))),sum(isnan(Angle_each_Controlline))];


 observed_freeview = [Frewview_Stim;Frewview_Control];


[p_freeview_prop,Q_freeview_prop] = chi2test(observed_freeview);

% Fixation 
Fixation_Stim = [length(Angle_each_Fix_line)-sum(isnan(Angle_each_Fix_line)),sum(isnan(Angle_each_Fix_line))];
Fixation_Control = [(length(Angle_each_Fix_Controlline)-sum(isnan(Angle_each_Fix_Controlline))),sum(isnan(Angle_each_Fix_Controlline))];


 observed_fixation = [Fixation_Stim;Fixation_Control];


[p_fixation_prop,Q_fixation_prop] = chi2test(observed_fixation);

% Gap Fixation
GapFixation_Stim = [length(Angle_each_GapFix_line)-sum(isnan(Angle_each_GapFix_line)),sum(isnan(Angle_each_GapFix_line))];
GapFixation_Control = [(length(Angle_each_GapFix_Controlline)-sum(isnan(Angle_each_GapFix_Controlline))),sum(isnan(Angle_each_GapFix_Controlline))];


 observed_gapfixation = [GapFixation_Stim;GapFixation_Control];


[p_gapfixation_prop,Q_gapfixation_prop] = chi2test(observed_gapfixation);


%SaccadeVector = [33.7500000000000	56.2500000000000	78.7500000000000	101.250000000000	123.750000000000	146.250000000000	168.750000000000	191.250000000000	213.750000000000	236.250000000000	258.750000000000	281.250000000000	303.750000000000	326.250000000000	348.750000000000	11.2500000000000];


%NumberOfAngle_Control_Sem = NaN*NumAngleControl';
%NumberOfAngle_Stim_Sem = NaN*NumAngleStim';
SeparationPoint=[22.5:22.5:360];
 
SaccadeVectorTuning=PoloarHistGroup(Angle_each_line_valid,ones(size(Angle_each_line_valid,2),1),SeparationPoint);
SaccadeVector=SaccadeVectorTuning.Center;
PropAngleStim=SaccadeVectorTuning.DataNum/sum(SaccadeVectorTuning.DataNum)*100;%Change into proportion of saccades
NumberOfAngle_Stim_Sem=NaN*ones(length(PropAngleStim),1);


SaccadeVectorTuning=PoloarHistGroup(Angle_each_Fixline_valid,ones(size(Angle_each_Fixline_valid,2),1),SeparationPoint);
SaccadeVector=SaccadeVectorTuning.Center;
PropAngleStim_Fix=SaccadeVectorTuning.DataNum/sum(SaccadeVectorTuning.DataNum)*100;%Change into proportion of saccades
NumberOfAngle_FixStim_Sem=NaN*ones(length(PropAngleStim_Fix),1);

SaccadeVectorTuning=PoloarHistGroup(Angle_each_GapFixline_valid,ones(size(Angle_each_GapFixline_valid,2),1),SeparationPoint);
SaccadeVector=SaccadeVectorTuning.Center;
PropAngleStim_GapFix=SaccadeVectorTuning.DataNum/sum(SaccadeVectorTuning.DataNum)*100;%Change into proportion of saccades
NumberOfAngle_GapFixStim_Sem=NaN*ones(length(PropAngleStim_GapFix),1);

%%Plotting
if ShowFigureFlag
    %Figure 1: Saccade distribution among aong the whole trial periods
    FigureStartNum=20;
    FigureIndex=1;
    figtitlestr{FigureIndex}='Distribution_Of_saccade_angleForThreeTasks';

    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
%

   % set(fig{FigureIndex}, 'PaperUnits', 'inches');
   % formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
%}
    %Control
   % pax = polaraxes;
   % polarplot(pax, [SaccadeVector/180*pi SaccadeVector(1)/180*pi],[NumAngleControl NumAngleControl(1)], '-r', 'LineWidth', 2);
   h = polarplot([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim PropAngleStim(1)],...
    'Color','r','LineWidth',1);
   
%h.MarkerSize = 0.5;

%stim
hold on

%s = polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim_Fix PropAngleStim_Fix(1)],[NumberOfAngle_Stim_Sem' NumberOfAngle_Stim_Sem(1)],...
 %   '-y','Color',[255,128,0]/255,'LineWidth',3);

s = polarplot([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim_Fix PropAngleStim_Fix(1)],...
    'Color',[255,128,0]/255,'LineWidth',1);

hold on

k = polarplot([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim_GapFix PropAngleStim_GapFix(1)],...
    'Color',[255,255,51]/255,'LineWidth',1);
%s.MarkerSize = 0.5;





set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',30,'FontWeight','Bold','LineWidth',1);
%set(gca,'FontSize',30,'LineWidth',0.5);


thetaticks([0:45:359]);
thetaticklabels([0:45:359]);

legend({'FreeView','Fix','GapFix'});




%thetaticks(sort(SaccadeVector));
%thetaticklabels(sort(SaccadeVector));

%title("Saccade Distribution",'FontSize',25,'FontWeight','Bold');


%rlim([0,max(max([NumAngleStim,NumAngleControl])*1.2)]);
%% Comparison of saccade amplitude
 FigureStartNum = FigureStartNum + 1;
 FigureIndex = FigureIndex + 1;
figtitlestr{FigureIndex}='SaccadeAmplitudeForThreeTasks';

   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1 1], 'nature');

 

b = bar([1,2,3],[Amp_free_mean,Amp_gapfix_mean,Amp_fix_mean]);
b.FaceColor='flat';
b.CData=[1,0,0;[255,128,0]/255;[255,255,51]/255];

hold on
errorbar([1,2,3],[Amp_free_mean,Amp_gapfix_mean,Amp_fix_mean],[Amp_free_sem,Amp_gapfix_sem,Amp_fix_sem],'.k');
xticklabels({'Free','GapFix','Fix'});
box off
ylabel('Saccade Amplitude(degree)');
set(gca,'LineWidth',0.5,'FontSize',5);

%% Comparison of saccade probability
 FigureStartNum = FigureStartNum + 1;
 FigureIndex = FigureIndex + 1;

figtitlestr{FigureIndex}='SaccadeProbabilityForThreeTasks';

   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});

    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [3 1.5], 'nature');


   subplot(1,3,1)
   b=bar([1,2],[Prop_Sac,Prop_Sac_Control]);
   b.FaceColor=['flat'];
   b.CData=[1,0,0;0,0,0];
   ylim([0,1]);
   box off;
   xticklabels({'Opto-Stim','Control'});
   ylabel('Saccade Probability');
   title('Free View');
   set(gca,'FontSize',5);

   subplot(1,3,2)
   b=bar([1,2],[Prop_Sac_GapFix,Prop_Sac_GapFix_Control]);
    ylim([0,1]);
    b.FaceColor=['flat'];
    b.CData=[[255,255,51]/255;0,0,0];
   ylim([0,1]);
   box off
   xticklabels({'Opto-Stim','Control'});
   title('Gap Fix');
   set(gca,'FontSize',5);

   subplot(1,3,3)
   b=bar([1,2],[Prop_Sac_Fix,Prop_Sac_Fix_Control]);
    ylim([0,1]);
b.FaceColor=['flat'];
    b.CData=[[255,128,0]/255;0,0,0];
   ylim([0,1]);
   box off
   xticklabels({'Opto-Stim','Control'});
   title('Fix');
   set(gca,'FontSize',5);

   


keyboard

end

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

%OutputFigureName='SC_VideoVeiw_Behavior';
OutputFigureName='FEF_FreeFixGapFix_Behavior';

if OutputFlag == 1
for jj=1:FigureIndex
        if ~isempty(fig{jj})
        

         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.pdf');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
end

    disp('Figures have been exported to the exported folder');
end
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

end

