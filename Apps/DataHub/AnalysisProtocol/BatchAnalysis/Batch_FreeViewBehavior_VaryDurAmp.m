function Batch_FreeViewBehavior_VaryDurAmp(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
FilesName_Each=Data.FileName;

UseAll = 1;%1: Use all the recorded sites; 0: only use significant sites 

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
     
     



     
     
     Task(i)=OutputData.FreeViewVaryDurAmpBehavior(1).TaskCode;
     TaskName{i}=OutputData.FreeViewVaryDurAmpBehavior(1).Task;
     
     
     Data =OutputData.FreeViewVaryDurAmpBehavior.DataStamp;
    


      Monkey = FilesName_Each{i}(1);

          if Monkey == 'R'
              MonkeyIndex(i) = 1;
          else
              MonkeyIndex(i) = 2;
          end


         VaryDur(i,:) = Data('Export_Dur');
        VaryPower(i,:) = Data('Export_Power');
         Prop_Sac_d_Control(i,:) = Data('Export_Prop_Sac_d_Control');
         Prop_Sac_d_Stim(i,:) = Data('Export_Prop_Sac_d_Stim');
        Prop_Sac_p_Control(i,:) = Data('Export_Prop_Sac_p_Control');
        Prop_Sac_p_Stim(i,:) = Data('Export_Prop_Sac_p_Stim')
        p_d(i,:) = Data('Export_p_d');
        p_p(i,:) = Data('Export_p_p');


   
     
      if isfield(OutputData,'FreeViewBehavior')

           Data =OutputData.FreeViewBehavior.DataStamp;
           p_prop(i) = Data('p_prop');

      end


     
     
     
end

 
Valid_dur = sum(isnan(p_d),2)==0;
Valid_p = sum(p_p<0.05,2)>0;

Valid_use_d = Valid_dur & Valid_p;


Prop_p_Control = Prop_Sac_p_Control(Valid_p,:);
Prop_p_Stim = Prop_Sac_p_Stim(Valid_p,:);

Prop_p_Control_Mean = mean(Prop_p_Control,1);
Prop_p_Control_Sem = std(Prop_p_Control,[],1)/sqrt(sum(Valid_p));

Prop_p_Stim_Mean = mean(Prop_p_Stim,1);
Prop_p_Stim_Sem = std(Prop_p_Stim,[],1)/sqrt(sum(Valid_p));


Prop_d_Control = Prop_Sac_d_Control(Valid_use_d ,:);
Prop_d_Stim = Prop_Sac_d_Stim(Valid_use_d ,:);

Prop_d_Control_Mean = mean(Prop_d_Control,1);
Prop_d_Control_Sem = std(Prop_d_Control,[],1)/sqrt(sum(Valid_use_d));

Prop_d_Stim_Mean = mean(Prop_d_Stim,1);
Prop_d_Stim_Sem = std(Prop_d_Stim,[],1)/sqrt(sum(Valid_use_d));

Power = [0.5,1.5,4.5];
Dur = [10,50,100,200];
%%Plotting
if ShowFigureFlag
    %Figure 1: Saccade distribution among aong the whole trial periods
    FigureStartNum=20;
    FigureIndex=1;
    figtitlestr{FigureIndex}='VaryAmp';

    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');


     plot(Power,Prop_p_Control','-k','Color',[160,160,160]/255,'LineWidth',0.5);
    hold on
    



    plot(Power,Prop_p_Stim','-r','Color',[255,153,153]/255,'LineWidth',0.5);
    hold on
    errorbar(Power,Prop_p_Control_Mean,Prop_p_Control_Sem,'-k','Color',[0,0,0],'LineWidth',1);
    errorbar(Power,Prop_p_Stim_Mean,Prop_p_Stim_Sem,'-r','Color',[1,0,0],'LineWidth',1);

    xticks(Power);
     xlim([Power(1)-0.5,Power(3)+0.5]);
box off
xlabel('Laser Power (mW)');
ylabel('Prop of Contralateral Saccade - Prop of Ipsilateral Saccade');
set(gca,'FontSize',5,'LineWidth',0.5);


 FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;
    figtitlestr{FigureIndex}='VaryDur';

    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[80,100, 600,600],'Name',figtitlestr{FigureIndex});


    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');

     plot(Dur,Prop_d_Control','-k','Color',[160,160,160]/255,'LineWidth',0.5);
    hold on
    



    plot(Dur,Prop_d_Stim','-r','Color',[255,153,153]/255,'LineWidth',0.5);
    hold on
    errorbar(Dur,Prop_d_Control_Mean,Prop_d_Control_Sem,'-k','Color',[0,0,0],'LineWidth',1);
    errorbar(Dur,Prop_d_Stim_Mean,Prop_d_Stim_Sem,'-r','Color',[1,0,0],'LineWidth',1);
    xticks(Dur);
    xlim([Dur(1)-10,Dur(4)+10]);
box off
xlabel('Duration of Opto-stim (ms)');
ylabel('Prop of Contralateral Saccade - Prop of Ipsilateral Saccade');
set(gca,'FontSize',5,'LineWidth',0.5);


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

%OutputFigureName='SC_VideoVeiw_Behavior';
OutputFigureName='FEF_VideoVeiwVaryAmpDur_Behavior';


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



end