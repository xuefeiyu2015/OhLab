function Batch_Optogenetics_Control(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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





     % Data =OutputData.FreeViewNeural.DataStamp;
     Data =OutputData.FreeViewNeuralRed.DataStamp;
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

 
     index = index+1;   
    end %End of for file name
     
end


%Plot controller
PlotOptoEffect = 1;


%Analysis controller 
AnalysisOS = 0;% For Video Free View


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
p_Cri50=p_Cri50(ScreenCell);

%maxnum_ms = max(cellfun(@numel,PSTH_Tgt_RF));

SM_mean = mean(SM);
SM_sem = std(SM)/length(SM);

[h,p_fromzero]=ttest2(SM,0);




mean_psth = mean(PSTH_norm,1);
sem_psth = std(PSTH_norm,[],1)/sqrt(size(PSTH_norm,1));






if ShowFigureFlag

    FigureIndex=1;
    FigureStartNum=45;

    %figtitlestr{FigureIndex}='PSTH_OptoModulation_Control';
    figtitlestr{FigureIndex}='PSTH_OptoModulation_ControlRed';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [2 2], 'nature');


     AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        %{
        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        %}
        shadedErrorBar(PSTH_Opto_Time,mean_psth,sem_psth,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);
        ylim([-0.5,1]);
       % legend(sprintf('N=%d', ModuTypeCount(1)));
        
       % title('Population');
       % xlabel('Time from Opto-stim onset (ms)')
        ylabel('Normalized Firing Rate(Hz)');
        

        set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
      %  xlabel('Time from Opto-stim Onset (ms)');
      %  ylabel('Normalized Responses');

        box off;

    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1;


   % figtitlestr{FigureIndex}='Histogram_OptoMod_Control';
   figtitlestr{FigureIndex}='Histogram_OptoMod_ControlRed';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [2 2], 'nature');

    bin_sep=-20:1:20;
    b = histogram(SM, bin_sep,'Facecolor','w');
     set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
    xlabel('Modulation Strength');
    ylabel('Number of neurons');
    box off


    






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

OutputFigureName='PreInjection_Optostim_FEF';

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