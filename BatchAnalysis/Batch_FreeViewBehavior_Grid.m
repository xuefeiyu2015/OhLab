function Batch_FreeViewBehavior_Grid(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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

%Load the grid file
gridLoc = importdata('/Users/yux8/WorkRelated_YXF/Program_Matlab_Local/DataAnalysis/Results/Grid_FEF_BothRA_VideoFreeView_behavior.txt');


for i=1:length(FilesName)
  clear OutputData;
  
     load(FilesName{i}{1});
     
     OutputDataTemp=OutputData;
     
     
     Task(i)=OutputData.FreeViewBehavior(1).TaskCode;
     TaskName{i}=OutputData.FreeViewBehavior(1).Task;
     
     
     Data =OutputData.FreeViewBehavior.DataStamp;
    


      Monkey = FilesName_Each{i}(1);

          if Monkey == 'R'
              MonkeyIndex(i) = 1;
          else
              MonkeyIndex(i) = 2;
          end


     PropClSac(i) = Data('PropClSac');
     p_prop(i) = Data('p_prop');

     NumberOfAngle_Control(i,:) = Data('NumberOfAngle_Control');
     NumberOfAngle_Stim(i,:) = Data('NumberOfAngle_Stim');

     SaccadeCountStim_Sum(i,:) = Data('SaccadeCountStim');
     SaccadeCountControl_Sum(i,:) = Data('SaccadeCountControl');

     theta_stim{i}=Data('theta_Stim');
     amp_stim{i}=Data('Amp_Stim');
     theta_control{i}=Data('theta_Control');
     amp_control{i}=Data('Amp_Control');

      PreferVectorNum_Stim(i)=Data('PrefVectorNum_Stim');
      PreferVectorNum_Control(i)=Data('PrefVectorNum_Control');
      PreferVector_AmpNum_Stim(i)=Data('PrefVector_AmpNum_Stim');

      AveAmp_Stim(i)=median(amp_stim{i});
     




     
     
     TimeSeq = Data('TimeForSaccadeRaster');

     Data =OutputData.FreeViewNeural.DataStamp;



     SM_p(i) = Data('StimModulation_p');

     %%
 
end

    

%Significant case 
Total = SM_p < 0.05 |p_prop < 0.05;
Sig = p_prop < 0.05;
Sig =Sig';
NumSig = sum(Sig);
NumTotal = sum(Total);
disp(sprintf('Total: %d',NumTotal));
disp(sprintf('SigNum: %d',NumSig));

MonkeyIndex=MonkeyIndex';

if UseAll
    PropCS_total = PropClSac;
    PropCS_Sig = PropClSac;
  %  Sig = p_prop<1;
else
    PropCS_total = PropClSac(Total);
    PropCS_Sig = PropClSac(Sig);

    MonkeyIndex=MonkeyIndex(Sig);


    theta_stim=cell2mat(theta_stim(Sig))+90;
    amp_stim=cell2mat(amp_stim(Sig));

    gridLoc=gridLoc(Sig,:);

    PreferVectorNum_Stim=PreferVectorNum_Stim(Sig);
    PreferVectorNum_Control=PreferVectorNum_Control(Sig);
   % PreferVector_AmpNum_Stim=PreferVector_AmpNum_Stim(Sig);
    AveAmp_Stim=AveAmp_Stim(Sig);

end

UniqueGridLoc_r = unique(gridLoc(MonkeyIndex==1,1:2),'rows');
NumGrid_r = size(UniqueGridLoc_r,1);

for i = 1:NumGrid_r
   sel_r = gridLoc(:,1)==UniqueGridLoc_r(i,1) & gridLoc(:,2)==UniqueGridLoc_r(i,2) & MonkeyIndex==1;
   PreferVec_Stim_r{i}=PreferVectorNum_Stim(sel_r);
   PreferVecAmp_Stim_r{i}=AveAmp_Stim(sel_r);
  % Amp_Ave_Stim_r(i)=nanmean(AveAmp_Stim(sel_r));
   Sig_r(i)=sum(Sig(sel_r))>0;
   Sig_r_all{i} = Sig(sel_r);


end

UniqueGridLoc_a = unique(gridLoc(MonkeyIndex==2,1:2),'rows');
NumGrid_a = size(UniqueGridLoc_a,1);

for i = 1:NumGrid_a
   
   sel_a = gridLoc(:,1)==UniqueGridLoc_a(i,1) & gridLoc(:,2)==UniqueGridLoc_a(i,2) & MonkeyIndex==2;
   PreferVec_Stim_a{i}=PreferVectorNum_Stim(sel_a);
   PreferVecAmp_Stim_a{i}=AveAmp_Stim(sel_a);
  % Amp_Ave_Stim_a(i)=nanmean(AveAmp_Stim(sel_a));
   Sig_a(i)=sum(Sig(sel_a))>0;
   Sig_a_all{i} = Sig(sel_a);

end


if ShowFigureFlag
%% For Robin first
 FigureStartNum=30;
 FigureIndex=1;
 figtitlestr{FigureIndex}='Behavior_Grid_Robin';

 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [4 4], 'nature');


index = 1;
for i =1:NumGrid_r
    subplot(ceil(NumGrid_r/4),4,index);
    NumVec = length(PreferVec_Stim_r{i});
    if NumVec > 0
        for j = 1:NumVec
            
               vec =  PreferVec_Stim_r{i}(j);
               amp =  PreferVecAmp_Stim_r{i}(j);

              if Sig_r_all{i}(j)>0
                h_pc = polarplot([vec vec]/180*pi,[0 amp],'r');
              else
                 h_pc = polarplot([vec vec]/180*pi,[0 0]); 
                 
              end
                hold on
               
            
            

        end
        set(gca,'ThetaDir' , 'counterclockwise');
            set(gca,'ThetaZeroLocation','top')
            rlim([0,30]);
            thetaticklabels([]);
            rticklabels([]);

        title(sprintf('%d,%d',UniqueGridLoc_r(i,1),UniqueGridLoc_r(i,2)));

    else
        axis off;


    end

    
   
    index=index+1;


end


%% Now it's Adams turn
 FigureStartNum=FigureStartNum+1;
 FigureIndex=FigureIndex+1;
 figtitlestr{FigureIndex}='Behavior_Grid_Adams';

 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 600,600],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [4 4], 'nature');


index = 1;
for i =1:NumGrid_a
    subplot(ceil(NumGrid_a/4),4,index);
    NumVec = length(PreferVec_Stim_a{i});
    if NumVec > 0
        for j = 1:NumVec
            
               vec =  PreferVec_Stim_a{i}(j);
               amp =  PreferVecAmp_Stim_a{i}(j);

              if Sig_a_all{i}(j)>0
                h_pc = polarplot([vec vec]/180*pi,[0 amp],'r');
              else
                 h_pc = polarplot([vec vec]/180*pi,[0 0]); 
                 
              end
                hold on
               
            
            

        end
        set(gca,'ThetaDir' , 'counterclockwise');
            set(gca,'ThetaZeroLocation','top')
            rlim([0,30]);
            thetaticklabels([]);
            rticklabels([]);

        title(sprintf('%d,%d',UniqueGridLoc_a(i,1),UniqueGridLoc_a(i,2)));

    else
        axis off;


    end

    
   
    index=index+1;


end


end % If showFigure on 


%keyboard



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
OutputFigureName='FEF_VideoVeiw_Behavior_Grid';

if OutputFlag
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

