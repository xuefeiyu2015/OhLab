function  ExploreSaccadeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analyse exploretary saccades and choices in the explore tasks;
%Xuefei Yu 04202020
%Unpack all the data base and task settings
ECode;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
DataPath= Data.Path;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel


RandLocTrials=[608:608+5];
FixLocTrials=[477,8,766,771,1001,4];

BestScenes=[1063,1057,1056];


BestScenes=[1054
        1055
        1056
        1057
        1060
        1061
        1062
        1063];
%}
%BestScenes=[];

%TrialForEyeTracePlot=[766,610];%First is the fix loc, second is the rand loc
TrialForEyeTracePlot=[1,10];
PlotExampleFlag=0;

SelectiveTargetNum=1;%Select the period searching for the first target

FileName=Data.FileName;

TargetOnOffset=0;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;
%{
%Criterias for selecting saccades
V_Threshold=30;
ContinueBin=5;
ISI_Threshold=20;
%}
%BlockSepTrial=unique(SceneLib(:,1));
%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');

BlockType=Data.BlockType(:,:);

if RecordDate<datestr('04,03,2021') && TaskCode==FORAGING_EXPLORE_TASK_ITL
    fix_tmp=BlockType(:,3)>0;
    BlockType(fix_tmp,2)=BlockType(fix_tmp,2)-10;
    
end
BlockType=BlockType(Selection,2);





%Correct early block type
%Before 03,18,2020, I only use 0 to indicate fix scene sequence and 1 to
%indicate random scene sequence, both with fixed location sequence
if RecordDate<datestr('03,17,2020')
    
   BlockType(BlockType==0)=10;
    BlockType(BlockType==1)=11;
   
    
end
%{
BlockType(BlockType==10)=11;
BlockType(BlockType==0)=1;
%}

%BlockType(BlockType(:,3)>0,3)=999;

[UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
NumberOfBlocks=max(reprojectorderindex);

%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

EyeBinWidth=Data.EyeBinWidth;

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);



%Recalibrate eyes according to the extra FP
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

ParaLib=Data.ParaLib;
ExtraFP_Info=ParaLib('ExtraFixPointID,Style,Size,Loc_X,Loc_Y');
ExtraFP_Loc=ExtraFP_Info(Selection,end-1:end);


%Get the information about the extra FP
%ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
ExtraEyeIn=FindOutTime(EventChannel,EventBin,EXTRAFIXSTARTCD);
ExtraEyeOff=FindOutTime(EventChannel,EventBin,EXTRAFIXOFF);

EyeInBin_X=SelectEyeInterval(EyeChannel_X_L,ExtraEyeIn,ExtraEyeOff);
EyeInBin_X=EyeInBin_X{1};
EyeInBin_Y=SelectEyeInterval(EyeChannel_Y_L,ExtraEyeIn,ExtraEyeOff);
EyeInBin_Y=EyeInBin_Y{1};
EyeDataForCal={EyeInBin_X,EyeInBin_Y};

EyeDataRaw={EyeChannel_X_L,EyeChannel_Y_L};
CalPlotCheck=0;%Examine the calibration results
CaledEye=EyeCalibration(EyeDataRaw,EyeDataForCal,EyeCalInfo_L,ExtraFP_Loc,CalPlotCheck);




EyeChannel_X_L=CaledEye.EyeTrace{1};
EyeChannel_Y_L=CaledEye.EyeTrace{2};



%Get scene and fractals information
SceneName=Data.SceneIndex(Selection,3);
TargetName=Data.TargetIndex(Selection,2:end,2);
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName),2);

TargetSize=Data.TargetSize(Selection,2:end);
TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);
TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);

ChoiceOutcome=Data.ChoiceOutcome(Selection,2:end);

%Select the trials at lease has one case of advance eye in

AdvanceInTargetNum=sum(EventChannel==ADVEYE_IN,2);
%AdvanceInTargetTrial=AdvanceInTargetNum>0;
AdvanceInTargetTrial=AdvanceInTargetNum==TargetNum;
strAdvanceIn=fprintf('There are %d Advance-Eye-In trials\n',sum(AdvanceInTargetTrial));

%Select Searching time
 FixOffTime=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);
 TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
 OutScreenTime=FindOutTime(EventChannel,EventTimeChannel,EYEOUTSCREEN);
 
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 TargetOnTime=TargetOnTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 
 
 
 FixOffTime=ScreenOutOverScreen(FixOffTime,OutScreenTime);
 
  %Remove the col
 NaNCol_Check=sum(isnan(TargetOnTime),1);

 TargetOnTime(:,NaNCol_Check==size(TargetOnTime,1))=[];

 
 TargetName=TargetName(AdvanceInTargetTrial,:,:);
 TargetMat=floor(~isnan(TargetName));
 TargetMat(TargetMat==0)=NaN;
 
 SceneName=SceneName(AdvanceInTargetTrial,:);
 
 EyeChannel_X_L=EyeChannel_X_L(AdvanceInTargetTrial);
 EyeChannel_Y_L=EyeChannel_Y_L(AdvanceInTargetTrial);
 
 TargetLoc=TargetLoc(AdvanceInTargetTrial,:,:);
 TargetPolarAngle=TargetPolarAngle(AdvanceInTargetTrial,:);
 TargetNum=TargetNum(AdvanceInTargetTrial);
 
 BlockType=BlockType(AdvanceInTargetTrial);
 
 
SearchingTime= TargetOnTime-FixOffTime;

ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FLocFScene','FLocRScene','RLocFScene','RLocRScene'};

%Analysis saccades using cluster fix

EyeChannel_reorganize=cellfun(@(x,y) [x;y],EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);

%Using clusterfix to analysis eyemovements
%It took a long time to form clusters.
%To save time, save the results to the data folder

path=DataPath;
cd(path);


SaccadeFileName=strcat(FileName(1:end-4),'_saccade','.mat');

if exist(SaccadeFileName)
    load(SaccadeFileName);
    disp('preanalysized saccades has been loaded to the workspace') ;
else
   disp('Analyzing the data using cluster fix, please wait.') 
   fixationstats  =ClusterFix(EyeChannel_reorganize,1/(EyeBinWidth*1000));
   save(SaccadeFileName,'fixationstats');
   disp('Saccades sucessfully exported'); 
   load(SaccadeFileName);
   disp('preanalysized saccades has been loaded to the workspace') ;

    
    
end

SaccadesTime=cellfun(@(x) x.saccadetimes,fixationstats ,'uniform',0);
Saccade_Value=cellfun(@(x) x.SaaccadeClusterValues,fixationstats ,'uniform',0);
Saccade_MaxAcceleration=cellfun(@(x) x(:,1),Saccade_Value ,'uniform',0);
Saccade_MaxV=cellfun(@(x) x(:,2),Saccade_Value ,'uniform',0);

Saccade_MeanDist=cellfun(@(x) x(:,3),Saccade_Value ,'uniform',0);
Saccade_MeanVelocity=cellfun(@(x) x(:,4),Saccade_Value ,'uniform',0);

Saccade_MeanAbsAngle=cellfun(@(x) x(:,5),Saccade_Value ,'uniform',0);
Saccade_MeanAngV=cellfun(@(x) x(:,6),Saccade_Value ,'uniform',0);

FixationTime=cellfun(@(x) x.fixationtimes,fixationstats ,'uniform',0);




%%%%%%%%%%%%%%%%%%%%%%%%First, analysis the saccade for each explortary period%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StartTimeBin=floor((FixOffTime+StartDelay)/EyeBinWidth);         
EndTimeBin=ceil((TargetOnTime+TargetOnOffset)/EyeBinWidth);

EyeChannel_X_Explore=SelectEyeInterval(EyeChannel_X_L,StartTimeBin,EndTimeBin);
EyeChannel_Y_Explore=SelectEyeInterval(EyeChannel_Y_L,StartTimeBin,EndTimeBin);


SearchingTime= TargetOnTime-FixOffTime;

%Saccade=SelectSaccadeFromCell(TrialIndex,EyeDataX,EyeDataY,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);



SearchingTime_Select=SearchingTime(:,SelectiveTargetNum);
EyeChannel_X_Explore_Select=EyeChannel_X_Explore{SelectiveTargetNum};
EyeChannel_Y_Explore_Select=EyeChannel_Y_Explore{SelectiveTargetNum};

StartTimeBin_Select=StartTimeBin(:,SelectiveTargetNum);
EndTimeBin_Select=EndTimeBin(:,SelectiveTargetNum);

%Measure how many saccades are found in each trial during the first explore
%period
SaccadesTime_Explore=SelectPeriod(SaccadesTime,[StartTimeBin_Select,EndTimeBin_Select],1);
FixTime_Explore=SelectPeriod(FixationTime,[StartTimeBin_Select,EndTimeBin_Select],0);



NumOfStep=cellfun(@(x) size(x,2),SaccadesTime_Explore);

FirstSaccadeStartTime=cellfun(@(x) x(1,1),SaccadesTime_Explore,'uniform',0);
FirstSaccadeEndTime=cellfun(@(x) x(2,1),SaccadesTime_Explore,'uniform',0);



ResponseLatency=(cell2mat(FirstSaccadeStartTime)'-StartTimeBin_Select)*EyeBinWidth;
ResponseLatency(ResponseLatency<0)=0;


%Find out all the fixation position during exploration period
 for j=1:length(EyeChannel_X_Explore)
     
     StartTimeBin_EachFractal=StartTimeBin(:,j);
     EndTimeBin_EachFractal=EndTimeBin(:,j);
     
    SaccadesTime_Explore_EachFractal=SelectPeriod(SaccadesTime,[StartTimeBin_EachFractal,EndTimeBin_EachFractal],1);
    FixTime_Explore_EachFractal=SelectPeriod(FixationTime,[StartTimeBin_EachFractal,EndTimeBin_EachFractal],0)';
    FixLoc_Explore_EachFractal{j}=cellfun(@(x,y,z) FindFixLoc(y,z,x),FixTime_Explore_EachFractal,EyeChannel_X_L,EyeChannel_Y_L,'uniform',0); 
    
    NumOfStepEachFractal(:,j)=cellfun(@(x) size(x,2),SaccadesTime_Explore_EachFractal);
     
 end

NumOfStepEachFractal=NumOfStepEachFractal.*TargetMat;
NumOfStepAllFractal=sum(NumOfStepEachFractal,2,'omitnan' );
%{
 %Meansure the number of steps as a function of the fractal order
    for j=1:length(EyeChannel_X_Explore)
         EyeChannel_X_Explore_Order=EyeChannel_X_Explore{j};
         EyeChannel_Y_Explore_Order=EyeChannel_Y_Explore{j};
         StartTimeBin_Order=StartTimeBin(:,j);
         EndTimeBin_Order=EndTimeBin(:,j);
         
         SaccadesTime_Explore_temp=SelectPeriod(SaccadesTime,[StartTimeBin_Order,EndTimeBin_Order]);
         FixTime_Explore_temp=SelectPeriod(FixationTime,[StartTimeBin_Order,EndTimeBin_Order]);
          NumOfStep_Order{j}=cellfun(@(x) size(x,2),SaccadesTime_Explore_temp);
        
    end
%}



for i=1:NumberOfBlocks
    
     CurrentBlock=UniqueBlockType(i);
     Order(i)=find(AssociateBlock==UniqueBlockType(i));
    %Select out the eye movements during the searching period for each
     %scene 
     SceneName_CurrBlock=SceneName(BlockType==CurrentBlock,:);
     UniqueScene_CurrBlock=unique(SceneName_CurrBlock);
     NumOfStep_currBlock=NumOfStep(BlockType==CurrentBlock);
     NumOfStep_currBlock_mean(i)=mean(NumOfStep_currBlock);
     NumOfStep_currBlock_sem(i)=std(NumOfStep_currBlock)/sqrt(sum(BlockType==CurrentBlock));
     
     NumOfStepAllFractal_currBlock=NumOfStepAllFractal(BlockType==CurrentBlock);
     NumOfAllFractal_currBlock=TargetNum(BlockType==CurrentBlock);
     
     NumOfStepAllFractal_currBlock_mean(i)=sum(NumOfStepAllFractal_currBlock.*NumOfAllFractal_currBlock)/sum(NumOfAllFractal_currBlock);
     NumOfStepAllFractal_currBlock_std=sqrt(sum((NumOfStepAllFractal_currBlock-NumOfStepAllFractal_currBlock_mean(i)).^2.*NumOfAllFractal_currBlock)/sum(NumOfAllFractal_currBlock));
     NumOfStepAllFractal_currBlock_sem(i)=NumOfStepAllFractal_currBlock_std./sqrt(sum(NumOfAllFractal_currBlock));
     
     ResponseLatency_currBlock_mean(i)=nanmean(ResponseLatency(BlockType==CurrentBlock));
     ResponseLatency_currBlock_sem(i)=nanstd(ResponseLatency(BlockType==CurrentBlock))/sqrt(sum(BlockType==CurrentBlock));
     
     ResponseLatency_Stats{i}=ResponseLatency(BlockType==CurrentBlock);
     NumOfStep_currBlock_Stats{i}=NumOfStep_currBlock;
     
     %Analysis the fixation location during exploration period
   %  FixLoc_Explore_EachFractal_currBlock_temp=cellfun(@(x) x(BlockType==CurrentBlock & ismember(SceneName,BestScenes)),FixLoc_Explore_EachFractal,'uniform',0);
    FixLoc_Explore_EachFractal_currBlock_temp=cellfun(@(x) x(BlockType==CurrentBlock),FixLoc_Explore_EachFractal,'uniform',0);
     
      

     temp=cellfun(@(x) cell2mat(cellfun(@(x2) cell2mat(x2), x,'uniform',0)'),FixLoc_Explore_EachFractal_currBlock_temp,'uniform',0);
 
     FixLoc_Explore_EachFractal_currBlock{i}=temp;
     FixLoc_Explore_AllFractal_currBlock{i}=cell2mat(temp);
     
     FixDistToCenter_EachFractal_currBlock{i}=cellfun(@(x) sqrt(sum(x.^2,1)),temp,'uniform',0);
     FixDistToCenter_AllFractal_currBlock{i}=cell2mat(cellfun(@(x) sqrt(sum(x.^2,1)),temp,'uniform',0));
%     NumOfStep_Order_currBlock=cellfun(@(x) x(BlockType==CurrentBlock), NumOfStep_Order,'uniform',0);
    
  %   NumOfStep_Order_currBlock_mean(i,:)=cellfun(@nanmean , NumOfStep_Order_currBlock,'uniform',0);
  %   NumOfStep_Order_currBlock_sem(i,:)=cellfun(@(x) nanstd(x)/sqrt(sum(BlockType==CurrentBlock)), NumOfStep_Order_currBlock,'uniform',0);
    %{
    for j=1:length(UniqueScene_CurrBlock)
        
        NumOfStep_SceneBlock(i,j)=mean(NumOfStep_currBlock(SceneName_CurrBlock==UniqueScene_CurrBlock(j)));
      
       
    end
     %}
   
     
    
   
end







%Example saccade and fixation distribution for a random loc condition and a
%fix loc condition

if PlotExampleFlag
    
    for i=1:length(TrialForEyeTracePlot)
        Eye_X{i}=EyeChannel_X_L{TrialForEyeTracePlot(i)};
        Eye_Y{i}=EyeChannel_Y_L{TrialForEyeTracePlot(i)};
        FixTime_Example{i}=FixationTime{TrialForEyeTracePlot(i)};
        
        FixTime_Explore_Example{i}=FixTime_Explore{TrialForEyeTracePlot(i)};
         Eye_X_explore{i}=EyeChannel_X_Explore_Select{TrialForEyeTracePlot(i)};
        Eye_Y_explore{i}=EyeChannel_Y_Explore_Select{TrialForEyeTracePlot(i)};
        
        SaccadeTime_Example{i}=SaccadesTime_Explore{TrialForEyeTracePlot(i)};
        EndTimeBin_explore_example(i)=EndTimeBin(TrialForEyeTracePlot(i));
        
        TargetLoc_example(i,:)=TargetLoc(TrialForEyeTracePlot(i),SelectiveTargetNum,:);
        
    end
end
%{
[p1_RL,~]=ranksum(ResponseLatency_Stats{1},ResponseLatency_Stats{2});
[p2_SS,~]=ranksum(NumOfStep_currBlock_Stats{1},NumOfStep_currBlock_Stats{2});
%}


     

if ShowFigureFlag
    
    titlestr{1}='BarForSaccadeSteps';

fig{1}=PrepareFigure(120,'w',[200,100, 800,600],'Name',titlestr{1});
subplot(2,1,1)
b=bar(NumOfStep_currBlock_mean);
%b=bar(NumOfStepAllFractal_currBlock_mean);

b.FaceColor = 'flat';
b.EdgeColor = 'flat';

b.CData = ColorForEachMode(Order,:);
hold on
box off
errorbar(NumOfStep_currBlock_mean,NumOfStep_currBlock_sem,'ok','LineWidth',3);
%errorbar(NumOfStepAllFractal_currBlock_mean,NumOfStepAllFractal_currBlock_sem,'ok','LineWidth',3);
ylabel('Average number of steps');
xticks=[1:NumberOfBlocks];
xticklabels(Blockstr(Order));
set(gca,'LineWidth',3,'FontSize',25);

subplot(2,1,2)
b=bar(ResponseLatency_currBlock_mean);
b.FaceColor = 'flat';
b.EdgeColor = 'flat';

b.CData = ColorForEachMode(Order,:);
hold on
box off
errorbar(ResponseLatency_currBlock_mean,ResponseLatency_currBlock_sem,'ok','LineWidth',3);
ylabel('Average first saccade latency');
xticks=[1:NumberOfBlocks];
xticklabels(Blockstr(Order));
set(gca,'LineWidth',3,'FontSize',25);





%Plot two example exploration saccades
if PlotExampleFlag
    
     titlestr{2}='ExampleSaccadesAndFixations';

     fig{2}=PrepareFigure(121,'w',[200,100, 800,300],'Name',titlestr{2});

     for i=1:length(TrialForEyeTracePlot)
         subplot(ceil(length(TrialForEyeTracePlot)/2),2,i);
         
          %Plot the region for scene 
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
          c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
         
         hold on
         
         %Overlay object location
         ColorFractal=[0.8500, 0.3250, 0.0980];
             
                r=7.5;
                t=cplot(r,TargetLoc_example(i,1),TargetLoc_example(i,2));
                t.FaceColor=ColorFractal;
                t.FaceAlpha=0.1;
                t.EdgeColor=ColorFractal;
                t.EdgeAlpha=0.1;
           
         
         %
         
         plot(Eye_X_explore{i},Eye_Y_explore{i},'-r','LineWidth',3);
         
         z=1:numel(Eye_X_explore{i});
         
         
         
        
          plot(Eye_X{i}(min(FixTime_Explore_Example{i}(1,:)):max(EndTimeBin_explore_example(i),max(FixTime_Explore_Example{i}(2,:)))),Eye_Y{i}(min(FixTime_Explore_Example{i}(1,:)):max(EndTimeBin_explore_example(i),max(FixTime_Explore_Example{i}(2,:)))),'-r','LineWidth',2);
         hold on
         hold on
         for j=1:size(FixTime_Explore_Example{i},2)
             temp=FixTime_Explore_Example{i};
             for jj=1:size(temp,2)
               plot(Eye_X{i}(temp(1,jj):temp(2,jj)),Eye_Y{i}(temp(1,jj):temp(2,jj)),'-b','LineWidth',3);
              % plot(Eye_X{i}(temp(2,jj)),Eye_Y{i}(temp(2,jj)),'-b','MarkerSize',10);
             end
             
             
             
         end
         
        
        % plot(Eye_X{i}(min(SaccadeTime_Example{i}(1,:)):max(SaccadeTime_Example{i}(2,:))),Eye_Y{i}(min(SaccadeTime_Example{i}(1,:)):max(SaccadeTime_Example{i}(2,:))),'-','LineWidth',3,'Color','#67726A');
        plot(Eye_X{i}(min(SaccadeTime_Example{i}(1,:)):max(SaccadeTime_Example{i}(2,:))),Eye_Y{i}(min(SaccadeTime_Example{i}(1,:)):max(SaccadeTime_Example{i}(2,:))),'-','LineWidth',3,'Color','r');
        hold on
         %Start Point
         plot(Eye_X{i}(min(SaccadeTime_Example{i}(1,:))),Eye_Y{i}(min(SaccadeTime_Example{i}(1,:))),'o','MarkerSize',10,'LineWidth',3,'Color','#11d646');
         %End Point
         plot(Eye_X{i}(max(SaccadeTime_Example{i}(2,:))),Eye_Y{i}(max(SaccadeTime_Example{i}(2,:))),'xg','MarkerSize',10,'LineWidth',3,'Color','#11d646');
         
         for j=1:size(SaccadeTime_Example{i},2)
             temp=SaccadeTime_Example{i};
             for jj=1:size(temp,2)
               plot(Eye_X{i}(temp(1,jj)),Eye_Y{i}(temp(1,jj)),'.b','MarkerSize',10);
               plot(Eye_X{i}(temp(2,jj)),Eye_Y{i}(temp(2,jj)),'.b','MarkerSize',10);
             end
             
             
             
         end
          for j=1:size(FixTime_Explore_Example{i},2)
             temp=FixTime_Explore_Example{i};
             for jj=1:size(temp,2)
               plot(Eye_X{i}(temp(1,jj):temp(2,jj)),Eye_Y{i}(temp(1,jj):temp(2,jj)),'.-b','MarkerSize',10);
               
             end
             
             
             
         end
         %}
         
         xlim([-30,30]);
         ylim([-30,30]);
         axis equal
         axis off
    
        
     end
    
    
 
else
    titlestr{2}=[];

     fig{2}=[];
     
     
end

%%%%%%%%%%%%%%%%% . Gaze location distribution%%%%%%%%%%%%%%%%%

 titlestr{3}='FixationDistribution';

fig{3}=PrepareFigure(122,'w',[200,100, 800,600],'Name',titlestr{3});  
for i=1:NumberOfBlocks
    
 
  subplot(ceil(NumberOfBlocks/2),2,i);
   %Plot the region for scene 
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
          c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
         
         hold on
         
         %Overlay object location
         ColorFractal=[0.8500, 0.3250, 0.0980];
             
                r=7.5;
                UniqueObjectLocWhole=unique(squeeze(TargetLoc(:,1,:)),'row');
               for ii=1:size(UniqueObjectLocWhole,1)
                t=cplot(r,UniqueObjectLocWhole(ii,1),UniqueObjectLocWhole(ii,2));
                t.FaceColor=ColorFractal;
                t.FaceAlpha=0.1;
                t.EdgeColor=ColorFractal;
                t.EdgeAlpha=0.1;
               end
           
         
         %Plot fixation location
         try
         currFixLocX=FixLoc_Explore_AllFractal_currBlock{i}(1,:);
         currFixLocY=FixLoc_Explore_AllFractal_currBlock{i}(2,:);
         catch
             currFixLocX=[];
             currFixLocY=[];
         end
         
         plot(currFixLocX,currFixLocY,'.b','MarkerSize',5);
         
         xlim([-30,30]);
         ylim([-30,30]);
         axis equal
         axis tight
         axis off
         
         title(Blockstr( Order(i)));
  
  
  
    
end

%%%%%%%%%%%%%%%%% . Gaze location distribution for each fractal order%%%%%%%%%%%%%%%%%

 titlestr{4}='FixationDistributionEachFractal';

fig{4}=PrepareFigure(123,'w',[200,100, 1000,1000],'Name',titlestr{4}); 
PPointer=1;
for i=1:NumberOfBlocks
   
    
    
    for j=1:length(EyeChannel_X_Explore)
    temp=FixLoc_Explore_EachFractal_currBlock{i}{j};
    if ~isempty(temp)
  subplot_tight(NumberOfBlocks,length(EyeChannel_X_Explore),PPointer,[0,0]);
 % figure;
   %Plot the region for scene 
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
          c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
         
         hold on
         
         %Overlay object location
         ColorFractal=[0.8500, 0.3250, 0.0980];
             
                r=7.5;
                UniqueObjectLocWhole=unique(squeeze(TargetLoc(:,1,:)),'row');
               for ii=1:size(UniqueObjectLocWhole,1)
                t=cplot(r,UniqueObjectLocWhole(ii,1),UniqueObjectLocWhole(ii,2));
                t.FaceColor=ColorFractal;
                t.FaceAlpha=0.1;
                t.EdgeColor=ColorFractal;
                t.EdgeAlpha=0.1;
               end
           
         
         %Plot fixation location
         currFixLocX=temp(1,:);
         currFixLocY=temp(2,:);
         
         plot(currFixLocX,currFixLocY,'.b','MarkerSize',5);
         
         xlim([-30,30]);
         ylim([-30,30]);
         axis equal
         axis off
         %ylabel(Blockstr( Order(i)))
         title(string(j));
         PPointer=PPointer+1;
    else
         PPointer=PPointer+1;
         
         
    end
    end
  
  
  
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%Distribution of gaze position distance for all
%%%%%%%%%%%%%%%%%%%%%%%%%%%fractals
titlestr{5}='FixationDistanceDistributionAllFractal';

fig{5}=PrepareFigure(124,'w',[200,100, 1000,400],'Name',titlestr{5}); 



for i=1:NumberOfBlocks
    
    dist=FixDistToCenter_AllFractal_currBlock{i};
    subplot(ceil(NumberOfBlocks/2),2,i);
    histogram(dist,'BinEdges',[0:0.5:25],'FaceColor',ColorForEachMode( Order(i),:));
    xlim([0,25]);
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    xlabel('Gaze Distance','FontSize',25,'FontWeight','Bold');
    ylabel('Number of gaze','FontSize',25,'FontWeight','Bold');
    title(Blockstr( Order(i)),'FontSize',25,'FontWeight','Bold');
    
   box off 
   % set(gca,'color','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%Distribution of gaze position distance for each
%%%%%%%%%%%%%%%%%%%%%%%%%%%fractals
titlestr{6}='FixationDistanceDistributionEachFractal';

fig{6}=PrepareFigure(125,'w',[200,100, 1000,600],'Name',titlestr{6}); 
PPointer=1;


for i=1:NumberOfBlocks
   
    
    
    for j=1:length(EyeChannel_X_Explore)
    dist=FixDistToCenter_EachFractal_currBlock{i}{j};
    if ~isempty(dist)
     % subplot_tight(NumberOfBlocks,length(EyeChannel_X_Explore),PPointer,[0.04,0.04]);
     subplot(NumberOfBlocks,length(EyeChannel_X_Explore),PPointer);
 
      
         
         %Plot fixation location distribution
         histogram(dist,'BinEdges',[0:0.5:25],'FaceColor',ColorForEachMode( Order(i),:),'EdgeColor',ColorForEachMode( Order(i),:));
    xlim([0,25]);
    xlabel('Gaze Distance');
    ylabel('Number of gaze');
   
         
         
         
         title(string(j));
         PPointer=PPointer+1;
    else
         PPointer=PPointer+1;
         
         
    end
    
    end
  
  
  
    
end


%%%%%%%%%%%Compare between fix loc and random loc in fixation exporation
titlestr{7}='FixationDistanceDistributionStrongComparison';

fig{7}=PrepareFigure(126,'w',[100,100, 1000,400],'Name',titlestr{7}); 


BlockType_II=BlockType;
BlockType_II(BlockType==0)=1;
BlockType_II(BlockType==10)=11;

BlockType_Original_Unique=unique(BlockType);

BlockType_II_unique=unique(BlockType_II);
BlockstrII={'RandomLoc','FixLoc'};

for i=1:length(BlockType_II_unique)
   Curr_BlockII=BlockType_II_unique(i);
   
    
   
    if Curr_BlockII==11
        dist=FixDistToCenter_EachFractal_currBlock(ismember(BlockType_Original_Unique,[10,11]));
        dist=cell2mat(cellfun(@(x) x{1},dist,'uniform',0));
        
    else
        dist=FixDistToCenter_EachFractal_currBlock(ismember(BlockType_Original_Unique,[0,1]));
        dist=cell2mat(cellfun(@(x) cell2mat(x),dist,'uniform',0));
    end
   
    %subplot_tight(ceil(NumberOfBlocks/2),2,i,[0.04,0.04]); 
    subplot(ceil(length(BlockType_II_unique)/2),2,i); 
    
         %Plot fixation location distribution
         g=histogram(dist,'BinEdges',[0:0.5:25],'FaceColor',ColorForEachMode(AssociateBlock==Curr_BlockII,:));
    xlim([0,25]);
    hold on
    yrange=ylim;
    plot([5,5],yrange,'--k','LineWidth',3);
    plot([7.5,7.5],yrange,'--k','LineWidth',3);
    plot([22.5,22.5],yrange,'--k','LineWidth',3);
    plot([15,15],yrange,'--k','LineWidth',3);
    
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    
    xlabel('Gaze Distance','FontSize',25,'FontWeight','Bold');
    ylabel('Number of gaze','FontSize',25,'FontWeight','Bold');
    title(BlockstrII{i},'FontSize',25,'FontWeight','Bold');
    
    
    box off;
   
    
    
    %{
    for j=1:length(EyeChannel_X_Explore)
    dist=FixDistToCenter_EachFractal_currBlock{i}{j};
    if ~isempty(dist)
     % subplot_tight(NumberOfBlocks,length(EyeChannel_X_Explore),PPointer,[0.04,0.04]);
     subplot(NumberOfBlocks,length(EyeChannel_X_Explore),PPointer);
 
      
         
         %Plot fixation location distribution
         histogram(dist,'BinEdges',[0:0.5:25],'FaceColor','k');
    xlim([0,25]);
    xlabel('Gaze Distance');
    ylabel('Number of gaze');
   
         
         
         
         title(string(j));
         PPointer=PPointer+1;
    else
         PPointer=PPointer+1;
         
         
    end
    end
  %}
  
  
    
end







end %End of the if plot flag

%{
trial=252;

Sac=SaccadesTime_Explore{trial};
fix=FixTime_Explore{trial};
figure
%plot(EyeChannel_X_Explore_Select{trial},EyeChannel_Y_Explore_Select{trial},'-r');
%hold on
plot(EyeChannel_X_L{trial},EyeChannel_Y_L{trial},'-r');

a=1;

%}

%{
for i=1:NumberOfBlocks
    
     CurrentBlock=UniqueBlockType(i);
    %Select out the eye movements during the searching period for each
     %scene 
     SceneName_CurrBlock=SceneName(BlockType==CurrentBlock,:);
     UniqueScene_CurrBlock=unique(SceneName_CurrBlock);
     
     SaccadeNumber_currBlock=arrayfun(@(x) x.NumOfSaccade,Saccades(BlockType==CurrentBlock));
     SaccadeStartTime=arrayfun(@(x) x.SaccadeStartTime,Saccades(BlockType==CurrentBlock),'uniform',0);
     FirstSaccadeStartTime=cellfun(@(x) x(1),SaccadeStartTime);
     
     SearchingTime_currBlock=SearchingTime(BlockType==CurrentBlock,:);
     SearchingTime_currBlock_mean(i)=nanmean(SearchingTime(BlockType==CurrentBlock));
     
    for j=1:length(UniqueScene_CurrBlock)
        
        FirstSaccadeStartTime_scene(i,j)=nanmean(FirstSaccadeStartTime(SceneName_CurrBlock==UniqueScene_CurrBlock(j)));
        FirstSaccadeStartTime_scene_sem(i,j)=nanstd(FirstSaccadeStartTime(SceneName_CurrBlock==UniqueScene_CurrBlock(j)))/sqrt(sum(SceneName_CurrBlock==UniqueScene_CurrBlock(j)));
        
        SearchingTime_scene(i,j)=nanmean(FirstSaccadeStartTime(SceneName_CurrBlock==UniqueScene_CurrBlock(j)));
        SearchingTime_scene_sem(i,j)=nanstd(FirstSaccadeStartTime(SceneName_CurrBlock==UniqueScene_CurrBlock(j)))/sqrt(sum(SceneName_CurrBlock==UniqueScene_CurrBlock(j)));
        
        SearchingTime_Scene_mean(i,j)=nanmean(SearchingTime_currBlock(SceneName_CurrBlock==UniqueScene_CurrBlock(j)));
        SearchingTime_Scene_sem(i,j)=nanstd(SearchingTime_currBlock(SceneName_CurrBlock==UniqueScene_CurrBlock(j)))/sqrt(sum(SceneName_CurrBlock==UniqueScene_CurrBlock(j)));
        
       
    end
     
      
     
     
    %{ 
    Saccades_currBlock=cellfun(@(x) x(BlockType==CurrentBlock),Saccades,'uniform',0);
    SaccadeNumber_currBlock=cell2mat(cellfun(@(x)  arrayfun(@(x2) x2.NumOfSaccade,x)',Saccades_currBlock,'uniform',0));
    SaccadeDuration_currBlock=cellfun(@(x)  arrayfun(@(x2) x2.SaccadeDuration,x,'uniform',0)',Saccades_currBlock,'uniform' ,0);
    SaccadeDurationWhole_currBlock=cell2mat(cellfun(@(x) cellfun(@(x2) sum(x2,'omitnan'),x),SaccadeDuration_currBlock,'uniform',0));
    
    SearchingTime_currBlock=SearchingTime(BlockType==CurrentBlock,SelectiveTargetNum);
    %}
    
   
end



figure
for i=1:NumberOfBlocks
subplot(2,1,1)

errorbar(FirstSaccadeStartTime_scene(i,:),FirstSaccadeStartTime_scene_sem(i,:),'linewidth',3);
hold on

subplot(2,1,2)
errorbar(SearchingTime_Scene_mean(i,:),SearchingTime_Scene_sem(i,:),'linewidth',3);
hold on
legend(string(UniqueBlockType(i)));


    
    
end

subplot(2,1,1)
legend(string(UniqueBlockType));

subplot(2,1,2)
legend(string(UniqueBlockType));

%}

%hold on
%errorbar(nanmean(SearchingTime_rep_Rand)',nanstd(SearchingTime_rep_Rand)','-sb','linewidth',2);
%%plot(SearchingTime_rep_Rand,'-ob','linewidth',2);


%legend({'RandLoc,FixSceneSeq','RandLoc,RandSceneSeq','FixLoc,FixSceneSeq','FixLoc,RandSceneSeq'}end