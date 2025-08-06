function  FreeViewEyeSaccadeSelectionIllustration(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);


ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;


if FileName(1)=='R'
    Monkey = 1;
else
    Monkey = 2;
end


Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%User handler
PlotExample = 1;%For plotting the raw for example, not needed for every session
%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);
Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);




EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

EyeTime = 0: EyeBinWidth:max(cellfun(@numel,EyeChannel_X_L))-EyeBinWidth;

EyeChannel_X_Line = cell2mat(reshape(EyeChannel_X_L,1,numel(EyeChannel_X_L)));
EyeChannel_Y_Line = cell2mat(reshape(EyeChannel_Y_L,1,numel(EyeChannel_Y_L)));


%Plot an eyetrace heatmap
screen_width = 40;
screen_height = 40;




[h, xedges, yedges] = histcounts2(EyeChannel_X_Line, EyeChannel_Y_Line, ...
    'XBinLimits', [-screen_width screen_width], ...
    'YBinLimits', [-screen_width screen_width], ...
    'NumBins', [80 80]);

%
% Apply Gaussian filter for smoothing
h_smooth_light = imgaussfilt(h, 2);  % Light smoothing
%h_smooth_heavy = imgaussfilt(h, 3);  % Heavy smoothing

imagesc(xedges(1:end-1), yedges(1:end-1), h_smooth_light');
colormap('jet');  % You can try other colormaps like 'hot', 'parula'
colorbar;
set(gca, 'YDir', 'normal');  % Ensure correct orientation
xlim([-40,40])
ylim([-40,40])
keyboard



%% Load the spike channel
SpikeChannelWhole=Data.SpikeChannel;
SpikeTimeWhole=Data.SpikeTimeData(:,:,Selection);
ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;

MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end

for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    
        SpikeTime=squeeze(SpikeTimeWhole(spk,:,:))';
        %{
    if SpikeChannel<4
        SpikeTime = SpikeTime;%Correct for the delay of online channel
    end
        %}
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);

    

%{
imagesc(xedges(1:end-1), yedges(1:end-1), h');
colormap('jet');  % You can try other colormaps like 'hot', 'parula'
colorbar;
set(gca, 'YDir', 'normal');  % Ensure correct orientation
xlim([-40,40])
ylim([-40,40])
keyboard

%Plot by time
EyeChannel_X_Line = cell2mat(reshape(EyeChannel_X_L,1,numel(EyeChannel_X_L)));
EyeChannel_Y_Line = cell2mat(reshape(EyeChannel_Y_L,1,numel(EyeChannel_Y_L)));

EyeTime_Line = 0: EyeBinWidth:length(EyeChannel_X_Line)-EyeBinWidth;


SelectT = EyeTime_Line>2000 & EyeTime_Line<3000;
Eye_X_SelectT = EyeChannel_X_Line(SelectT);
Eye_Y_SelectT = EyeChannel_Y_Line(SelectT);

Saccades=SelectSaccade(Eye_X_SelectT ,Eye_Y_SelectT ,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

%SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=cell2mat(arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)');
SaccadeAmplitude=cell2mat(arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)');

SaccadeStartTime_ori=cell2mat(arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)');
SaccadeEndTime_ori=cell2mat(arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)');

SaccadeStartPoint=cell2mat(arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)');
SaccadeEndPoint=cell2mat(arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)');
Sac_V = cell2mat(arrayfun(@(x)  x.PeakV,Saccades,'UniformOutput' ,0)');

SelectSac = SaccadeAmplitude> 5 & SaccadeAmplitude< 40 & Sac_V<1000;
Num = sum(SelectSac);
SaccadeAngle_Sel =SaccadeAngle(SelectSac);
SaccadeAmp_Sel =SaccadeAmplitude(SelectSac);
Sac_StartT = SaccadeStartTime_ori(SelectSac,:)+1;
Sac_EndT = SaccadeEndTime_ori(SelectSac,:)+1;
Sac_StartPoint = SaccadeStartPoint(SelectSac,:);
Sac_EndPoint = SaccadeEndPoint(SelectSac,:);

figure
%set(gca,'Color',[0.3,0.3,0.3]);
plot(Eye_X_SelectT,Eye_Y_SelectT,'-k','LineWidth',3);
hold on
for i = 1:Num
    
    plot(Eye_X_SelectT(Sac_StartT(i):Sac_EndT(i)),Eye_Y_SelectT(Sac_StartT(i):Sac_EndT(i)),'-r','LineWidth',3);
    hold on
    
end

xlim([-40,40]);
ylim([-40,40]);







keyboard



%Plot by trial

EyeChannel_X_Task=cell2mat(cellfun(@(x) x(1:2000),EyeChannel_X_L,'UniformOutput',false));
EyeChannel_Y_Task=cell2mat(cellfun(@(x) x(1:2000),EyeChannel_Y_L,'UniformOutput',false));





Trial =45;
Select = EyeTime>0 & EyeTime<2000;
Eye_X_Select = EyeChannel_X_Task(Trial,Select);
Eye_Y_Select = EyeChannel_Y_Task(Trial,Select);




Saccades=SelectSaccade(Eye_X_Select ,Eye_Y_Select ,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

%SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=cell2mat(arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)');
SaccadeAmplitude=cell2mat(arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)');

SaccadeStartTime_ori=cell2mat(arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)');
SaccadeEndTime_ori=cell2mat(arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)');

SaccadeStartPoint=cell2mat(arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)');
SaccadeEndPoint=cell2mat(arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)');
Sac_V = cell2mat(arrayfun(@(x)  x.PeakV,Saccades,'UniformOutput' ,0)');

SelectSac = SaccadeAmplitude> 5 & SaccadeAmplitude< 40;
Num = sum(SelectSac);
SaccadeAngle_Sel =SaccadeAngle(SelectSac);
SaccadeAmp_Sel =SaccadeAmplitude(SelectSac);
Sac_StartT = SaccadeStartTime_ori(SelectSac,:);
Sac_EndT = SaccadeEndTime_ori(SelectSac,:);
Sac_StartPoint = SaccadeStartPoint(SelectSac,:);
Sac_EndPoint = SaccadeEndPoint(SelectSac,:);


figure
%set(gca,'Color',[0.3,0.3,0.3]);
plot(Eye_X_Select,Eye_Y_Select,'-k','LineWidth',3);
hold on
for i = 1:Num
    
    plot(Eye_X_Select(Sac_StartT(i):Sac_EndT(i)),Eye_Y_Select(Sac_StartT(i):Sac_EndT(i)),'-r','LineWidth',3);
    hold on
end

xlim([-40,40]);
ylim([-40,40]);


%}
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
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end));%,OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
end
cd(OutPath);
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


%Output file name
OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
OutputFigureName=sprintf('%s_%s_N%s_C%s_%CID_%NID',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),string(ChannelID),string(ClusterID));





ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'FreeViewRawSaccade')
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

OutputData_New.FreeViewRawSaccade.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data


NameStr={'h','xedges', 'yedges'
    };

DataLib={h,xedges, yedges



};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.FreeViewRawSaccade.Task=TaskType;
OutputData_New.FreeViewRawSaccade.TaskCode=TaskCode;
OutputData_New.FreeViewRawSaccade.DataStamp=DataStamp;

OutputData_New.FreeViewRawSaccade.StartTrial =StartTrial;
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
Task_Old=OutputData.FreeViewRawSaccade.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
        %{
     elseif ismember(FileName,{'Adams062124VIDEOVIEW1'})
         %May save different task segment
         if StartTrial~=OutputData.FreeViewNeural.StartTrial

             if StartTrial > 160
                 OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(2),string(ChannelNum));
                 OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             elseif StartTrial > 360
                  OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(3),string(ChannelNum));
                  OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             end


         end
     elseif ismember(FileName,{'Robin030524VIDEOVIEW1'})
        % if StartTrial~=OutputData.FreeViewNeural.StartTrial

             if StartTrial > 80 && StartTrial <= 160
                 OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(2),string(ChannelNum));
                 OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;
 
             elseif StartTrial > 240 && StartTrial <= 160
                  OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(3),string(ChannelNum));
                  OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             end


       %  end
        %}     
     else
         %If the same task, replace the old one with the new one
         %OutputData=OutputData_New;
         OutputData.FreeViewRawSaccade=OutputData_New.FreeViewRawSaccade;
         
    
    
     end
else
    %Output the current file
   OutputData.FreeViewRawSaccade=OutputData_New.FreeViewRawSaccade; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  


end %Loop of spike
end

function organized=ReorganizeEye(data);
organized=[];
for i=1:length(data)
    data_curr=data{i};
    maxnumel=max(cellfun(@numel,data_curr));
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    
end

 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');

end

function ObjectIndex=ReproduceFromEvent(event,code)
for i=1:size(event,1)
    
    ObjectIndex{i}=event(i,ismember(event(i,:),code));
    
  
    
end
numeach=cellfun(@numel,ObjectIndex);
nummax=max(numeach);
if nummax==1
    ObjectIndex=cell2mat(ObjectIndex);

    
end


end