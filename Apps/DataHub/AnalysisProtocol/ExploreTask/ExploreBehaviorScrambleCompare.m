function ExploreBehaviorScrambleCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Unpack all the data base and task settings
ECode;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

TargetOnOffset=0;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;


%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);

BlockType=BlockType(Selection,2);

%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');
%To get the monkey name
RemainFiles=erase(Data.FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);

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

if RecordDate<datestr('10,28,2020') && strcmp(MonkeyName,'Robin')
try
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

catch
    disp('Skip offline eye calibration');
end

    
end


%Get scene and fractals information
SceneName=Data.SceneIndex(Selection,3);


TargetName=Data.TargetIndex(Selection,:,2);
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName(:,2:end)),2);
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
 %TargetOnTime=FixOffTime+300;
 AdvanceEyeInTime=FindOutTime(EventChannel,EventTimeChannel,ADVEYE_IN );
 
 TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
 OutScreenTime=FindOutTime(EventChannel,EventTimeChannel,EYEOUTSCREEN);
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 
 FixOffTime=ScreenOutOverScreen(FixOffTime,OutScreenTime);
 
 
 TargetOnTime=TargetOnTime(AdvanceInTargetTrial,:);
 %Remove the col
 NaNCol_Check=sum(isnan(TargetOnTime),1);

 TargetOnTime(:,NaNCol_Check==size(TargetOnTime,1))=[];
 
 SceneName=SceneName(AdvanceInTargetTrial,:);
 
 EyeChannel_X_L=EyeChannel_X_L(AdvanceInTargetTrial);
 EyeChannel_Y_L=EyeChannel_Y_L(AdvanceInTargetTrial);
 
 TargetLoc=TargetLoc(AdvanceInTargetTrial,:,:);
 TargetPolarAngle=TargetPolarAngle(AdvanceInTargetTrial,:);
 
 TargetNum=TargetNum(AdvanceInTargetTrial);
 BlockType=BlockType(AdvanceInTargetTrial);
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime)/1000;%Change into s for presentation

ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};


%Color Definition for each scenes
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};

%Setup scene library
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));
ImageFolerName='VisualStimulus';
%ImagePath=strcat(BasicDirectory,ImageFolerName,BasicDirectory(end));

if ShowFigureFlag
%Set up figure for eye position distribution
for i=1:NumberOfBlocks
if ishandle(101+i) 
    close(101+i); end
figure(101+i);set(101+i,'color','w');
Order=find(AssociateBlock==UniqueBlockType(i));
set(101+i,'Position', [100+i*10,100+i*5 1000,500], 'Name', strcat('Exploration for ',Blockstr{Order}));


end

if ishandle(106) 
    close(106); end
figure(106);set(106,'color','w');
set(106,'Position', [100+(i+1)*10,100+(i+1)*5 1000,1000], 'Name', 'Eye location Aligned for Each Mode');
end

ColorFractal=[0.8500, 0.3250, 0.0980];
Bin=0.5;
SepPoint=6; 
%Edges=[0:Bin:max(SearchingTime/1000)];
Edges=[0:Bin:SepPoint,max(max(max(SearchingTime)),SepPoint)+1];


%BlockTypeLine=[];
%FractalOrderLine=[];
for i=1:NumberOfBlocks
%for i=3:3
    CurrentBlock=UniqueBlockType(i);
    BlockOrder(i)=find(AssociateBlock==CurrentBlock);
    
    
    %Measure the average searching time for the current block
    SearchingTime_CurrBlock=SearchingTime(BlockType==CurrentBlock,:);
    SearchingTime_CurrBlock_Line{i}=reshape(SearchingTime_CurrBlock,1,numel(SearchingTime_CurrBlock));
    SearchingTime_CurrBlock_Count(i,:)= histcounts(SearchingTime_CurrBlock_Line{i},Edges);
    SearchingTime_CurrBlock_For_Stats{i}=SearchingTime_CurrBlock;
    
   % NumberOfTrials_CurrBlock(i)=sum(BlockType==CurrentBlock);
    
   % BlockTypeLine=[BlockTypeLine;repmat(CurrentBlock,NumberOfTrials_CurrBlock(i),1)];
   % FractalOrderLine=[FractalOrderLine;
    
    SearchingTime_MeanAll(i)=nanmean(SearchingTime_CurrBlock_Line{i});
    SearchingTime_SemAll(i)=nanstd(SearchingTime_CurrBlock_Line{i})./sqrt(sum(~isnan(SearchingTime_CurrBlock_Line{i})));
    
   
   
    SearchingTime_Mid{i}=median(SearchingTime_CurrBlock,1,'omitnan');
    CI_temp=ConfidenceInterval(SearchingTime_CurrBlock,90);
    SearchingTime_CI_Low{i}=CI_temp(1,:);
    SearchingTime_CI_High{i}=CI_temp(2,:);
    
    
    SearchingTime_Mean{i}=nanmean(SearchingTime_CurrBlock,1);
     SearchingTime_SEM_temp=nanstd(SearchingTime_CurrBlock,1)./sqrt(sum(~isnan(SearchingTime_CurrBlock),1));
     
     
     
     NumOfData_temp=sum(~isnan(SearchingTime_CurrBlock));
     SearchingTime_SEM_temp(NumOfData_temp<3)=NaN;
     
     SearchingTime_SEM{i}=SearchingTime_SEM_temp;
     
     
     
     
     
     NumOfData{i}=NumOfData_temp;
     NumberOfTrials_CurrBlock(i)=sum(NumOfData_temp);
     
     BlockIndex{i}=1:max(TargetNum,[],'all');
    
    
     
     
     
     
     
     
     %Select out the eye movements during the searching period for each
     %scene 
     SceneName_CurrBlock=SceneName(BlockType==CurrentBlock,:);
     UniqueScene_CurrBlock=unique(SceneName_CurrBlock);
     EyeChannel_X_CurrBlock=EyeChannel_X_L(BlockType==CurrentBlock);
     EyeChannel_Y_CurrBlock=EyeChannel_Y_L(BlockType==CurrentBlock);
     
     %Select out the object location
     TargetLoc_CurrBlock=TargetLoc(BlockType==CurrentBlock,:,:);
     TargetPolarAngle_CurrBlock=TargetPolarAngle(BlockType==CurrentBlock,:,:);
     
     
     
     FixOffTime_CurrBlock=FixOffTime(BlockType==CurrentBlock,:);
     TargetOnTime_CurrBlock=TargetOnTime(BlockType==CurrentBlock,:);
     
     
     %Separate Scenes into Scramble Scene and Normal Scenes
     ScrambleIndex_CurrBlock=SceneName_CurrBlock<0;
     NormalSceneIndex_CurrBlock=~ScrambleIndex_CurrBlock; 
     
     SearchingTimeScramble_CurrBlock=SearchingTime_CurrBlock(ScrambleIndex_CurrBlock,:);
     SearchingTimeNormal_CurrBlock=SearchingTime_CurrBlock( NormalSceneIndex_CurrBlock,:);
     
     SearchingTimeScramble_Mean{i}=nanmean(SearchingTimeScramble_CurrBlock,1);
     SearchingTimeScramble_SEM_temp=nanstd(SearchingTimeScramble_CurrBlock,1)./sqrt(sum(~isnan(SearchingTimeScramble_CurrBlock),1));
     

     NumOfData_temp=sum(~isnan(SearchingTimeScramble_CurrBlock));
     SearchingTimeScramble_SEM_temp(NumOfData_temp<3)=NaN;
 
     
     SearchingTimeScramble_SEM{i}=SearchingTimeScramble_SEM_temp;
     
   
     
     
      SearchingTimeNormal_Mean{i}=nanmean(SearchingTimeNormal_CurrBlock,1);
     SearchingTimeNormal_SEM_temp=nanstd(SearchingTimeNormal_CurrBlock,1)./sqrt(sum(~isnan(SearchingTimeNormal_CurrBlock),1));
     

     NumOfData_temp=sum(~isnan(SearchingTimeNormal_CurrBlock));
     SearchingTimeNormal_SEM_temp(NumOfData_temp<3)=NaN;
     
     SearchingTimeNormal_SEM{i}=SearchingTimeNormal_SEM_temp;
     
       %Use Matrix Representation
     %Only Fractal1 ,2 >2
     for iii=1:3
         if iii<3
             STScramble_Mat_Mean(i,iii)=nanmean(SearchingTimeScramble_CurrBlock(:,iii),'all');
             STScramble_Mat_SEM(i,iii)=nanstd(SearchingTimeScramble_CurrBlock(:,iii),[],'all')/sqrt(sum(~isnan(SearchingTimeScramble_CurrBlock(:,iii)),'all'));
             STScramble_Mat_Raw{i,iii}=SearchingTimeScramble_CurrBlock(:,iii);
             
             STNormal_Mat_Mean(i,iii)=nanmean(SearchingTimeNormal_CurrBlock(:,iii),'all');
             STNormal_Mat_SEM(i,iii)=nanstd(SearchingTimeNormal_CurrBlock(:,iii),[],'all')/sqrt(sum(~isnan(SearchingTimeNormal_CurrBlock(:,iii)),'all'));
             STNormal_Mat_Raw{i,iii}=SearchingTimeNormal_CurrBlock(:,iii);
             
             
         
         else
             STScramble_Mat_Mean(i,iii)=nanmean(SearchingTimeScramble_CurrBlock(:,iii:end),'all');
             STScramble_Mat_SEM(i,iii)=nanstd(SearchingTimeScramble_CurrBlock(:,iii:end),[],'all')/sqrt(sum(~isnan(SearchingTimeScramble_CurrBlock(:,iii:end)),'all'));
             STScramble_Mat_Raw{i,iii}=reshape(SearchingTimeScramble_CurrBlock(:,iii:end),numel(SearchingTimeScramble_CurrBlock(:,iii:end)),1);
             
              STNormal_Mat_Mean(i,iii)=nanmean(SearchingTimeNormal_CurrBlock(:,iii:end),'all');
             STNormal_Mat_SEM(i,iii)=nanstd(SearchingTimeNormal_CurrBlock(:,iii:end),[],'all')/sqrt(sum(~isnan(SearchingTimeNormal_CurrBlock(:,iii:end)),'all'));
             STNormal_Mat_Raw{i,iii}=reshape( SearchingTimeNormal_CurrBlock(:,iii:end),numel( SearchingTimeNormal_CurrBlock(:,iii:end)),1);
             
            
             
         end
         
     end
     
     
     
     
     
     
     for j=1:length(UniqueScene_CurrBlock)
     % for j=3:3
         %Select the eye trace for each scene
         Select_Scene=SceneName_CurrBlock==UniqueScene_CurrBlock(j);
         StartTimeBin=floor((FixOffTime_CurrBlock(Select_Scene,:)+StartDelay)/EyeBinWidth);         
         EndTimeBin=ceil((TargetOnTime_CurrBlock(Select_Scene,:)+TargetOnOffset)/EyeBinWidth);
         
         EyeChannel_X_CurrScene=EyeChannel_X_CurrBlock(Select_Scene);
         EyeChannel_Y_CurrScene=EyeChannel_Y_CurrBlock(Select_Scene);
         
         EyeChannel_X_Explore{j}=SelectEyeInterval(EyeChannel_X_CurrScene,StartTimeBin,EndTimeBin);
         EyeChannel_Y_Explore{j}=SelectEyeInterval(EyeChannel_Y_CurrScene,StartTimeBin,EndTimeBin);
         
         FractalOrder=1;
         
         EyeChannel_X_First_Explore=EyeChannel_X_Explore{j}{FractalOrder};
         EyeChannel_Y_First_Explore=EyeChannel_Y_Explore{j}{FractalOrder};
         
         %Select out current target location
         FirstTargetLoc_CurrScene=squeeze(TargetLoc_CurrBlock(Select_Scene,FractalOrder,:));
         TargetPolarAngle_CurrScene=squeeze(TargetPolarAngle_CurrBlock(Select_Scene,FractalOrder,:));
         if size(FirstTargetLoc_CurrScene,2)==1
             FirstTargetLoc_CurrScene=FirstTargetLoc_CurrScene';
         end
         
         
         
        % [xi,yi,pdfxy]=EyeDistDensity(EyeChannel_X_First_Explore,EyeChannel_Y_First_Explore);
         if ShowFigureFlag
         %Plot the first exploration eye trace for each scene
         figure(101+i);
          %Overlap images 
         CurrSceneName=UniqueScene_CurrBlock(j);
         SceneNameStr=sprintf('Scene%d.gif',CurrSceneName);
         
        % scene=imread(strcat(ImagePath,SceneNameStr));
        % figure
       %   axes('pos',[0 0 1 1])
        % imshow(scene);

         subplot(ceil(length(UniqueScene_CurrBlock)/4),4,j);
         
         %Plot the region for scene 
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
        % x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),'-');
         c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
     
         hold on
         f=cellfun(@plot,EyeChannel_X_First_Explore,EyeChannel_Y_First_Explore,'UniformOutput',0); 
         
          Eye_X_End=NaN*ones(length(EyeChannel_X_First_Explore),1);
          Eye_Y_End=Eye_X_End;
          %Empty check 
          empty_idx= cellfun(@(x) isempty(x),EyeChannel_X_First_Explore);

         
         Eye_X_End(~empty_idx)=cellfun(@(x) x(end),EyeChannel_X_First_Explore(~empty_idx));
         Eye_Y_End(~empty_idx)=cellfun(@(x) x(end),EyeChannel_Y_First_Explore(~empty_idx));
          
         plot(Eye_X_End,Eye_Y_End,'^', 'MarkerFaceColor',Color{j},'MarkerSize',5,'MarkerEdgeColor',Color{j})
        
         
         for ii=1:length(f)
        % f(ii).Marker='o';
            if ~isempty(f{ii})
             f{ii}.LineStyle='-';
             f{ii}.Color=Color{j};
             f{ii}.Color(4)=0.4;
             f{ii}.LineWidth=2;
             f{ii}.MarkerSize=1;
            end
         end
         
         %Overlay object location
       
         FirstTargetLoc_CurrSceneUnique=uniquetol(FirstTargetLoc_CurrScene,'ByRows',true);
         
         FirstTargetLoc_CurrScene_X=FirstTargetLoc_CurrSceneUnique(:,1);
         FirstTargetLoc_CurrScene_Y=FirstTargetLoc_CurrSceneUnique(:,2);
         
         
         
         for iii=1:length(FirstTargetLoc_CurrScene_X)
             if ~isnan(FirstTargetLoc_CurrScene_X(iii))
             
                r=7.5;
                t=cplot(r,FirstTargetLoc_CurrScene_X(iii),FirstTargetLoc_CurrScene_Y(iii));
                t.FaceColor=ColorFractal;
                t.FaceAlpha=0.1;
                t.EdgeColor=ColorFractal;
                t.EdgeAlpha=0.1;
             end
             
         end
         
         
         
         xlim([-30,30])
         ylim([-30,30])
        title(SceneNameStr(1:end-4),'FontWeight','Bold','FontSize',15);
        axis equal;
        axis off;
         end
       %Now shift the location of the target into a fixed location to overlay all eye trace 
       %Shift it anticlockwise
       %Set upward as zero
       
       temp_AfterTransformS_X=[];
       temp_AfterTransformS_Y=[];
       temp_AfterTransformN_X=[];
       temp_AfterTransformN_Y=[];
       
       temp_AfterTransform_X=[];
       temp_AfterTransform_Y=[];
           
       %Just Whether current scene is a scramble scene or not
       %{
      if UniqueScene_CurrBlock(j)<0
          %Scramble Scene
           temp_AfterTransformS_X=[];
           temp_AfterTransformS_Y=[];
          
      end
      %}
       
       for iii=1:length(TargetPolarAngle_CurrScene)
          if ~isnan(TargetPolarAngle_CurrScene(iii))
           TransformMat=[cos(TargetPolarAngle_CurrScene(iii)*pi/180),-sin(TargetPolarAngle_CurrScene(iii)*pi/180);sin(TargetPolarAngle_CurrScene(iii)*pi/180),cos(TargetPolarAngle_CurrScene(iii)*pi/180)];
       
           EyeChannel_X_First_Explore_Before=EyeChannel_X_First_Explore{iii};
           EyeChannel_Y_First_Explore_Before=EyeChannel_Y_First_Explore{iii};
           
           if isempty(EyeChannel_X_First_Explore_Before)
               temp_AfterTransform_X=[NaN];
               temp_AfterTransform_Y=[NaN];
               
               if UniqueScene_CurrBlock(j)<0
                   temp_AfterTransformS_X=[NaN];
                   temp_AfterTransformS_Y=[NaN];
                   
                   
               else
                   temp_AfterTransformN_X=[NaN];
                   temp_AfterTransformN_Y=[NaN];
                   
                   
               end
               
               
               
           else
                temp=[EyeChannel_X_First_Explore_Before',EyeChannel_Y_First_Explore_Before']*TransformMat;
                
                 temp_AfterTransform_X=[temp_AfterTransform_X,temp(:,1)',NaN];
                temp_AfterTransform_Y=[temp_AfterTransform_Y,temp(:,2)',NaN];
                
                if UniqueScene_CurrBlock(j)<0
                   temp_AfterTransformS_X=[temp_AfterTransformS_X,temp(:,1)',NaN];
                   temp_AfterTransformS_Y=[temp_AfterTransformS_Y,temp(:,2)',NaN];
                   
                   
               else
                   temp_AfterTransformN_X=[temp_AfterTransformN_X,temp(:,1)',NaN];
                   temp_AfterTransformN_Y=[temp_AfterTransformN_Y,temp(:,2)',NaN];
                   
               end
               
           
           end
       
    
          end
       
       end
       
       
       EyeChannel_X_Transform{j}=temp_AfterTransform_X;
       EyeChannel_Y_Transform{j}=temp_AfterTransform_Y;
       
        EyeChannelS_X_Transform{j}=temp_AfterTransformS_X;
       EyeChannelS_Y_Transform{j}=temp_AfterTransformS_Y;
       
       EyeChannelN_X_Transform{j}=temp_AfterTransformN_X;
       EyeChannelN_Y_Transform{j}=temp_AfterTransformN_Y;
       
       
     %{  
    %  figure(106);
    % subplot(2,4,j);

    % p=cellfun(@plot,EyeChannel_X_Transform,EyeChannel_Y_Transform)
    p=plot(temp_AfterTransform_X,temp_AfterTransform_Y)
     p.LineStyle='-';
         p.Color=Color{j};
         p.Color(4)=0.25;
         p.LineWidth=2;
         p.MarkerSize=1;
    
     xlim([-40,40])
     ylim([-40,40])
        
     %}   
        %  set(gca,'FontWeight','Bold');
        %set(gca,'FontSize',15);
        %set(gca,'LineWidth',2)
         
         
       
         
          %Plot the first exploration eye position density for each scene
       %  figure(201+i);
         
        % subplot(2,4,j);
        %{
        figure
         contourf(xi,yi,pdfxy);
         xlim([-40,40])
         ylim([-40,40])
         %}
  
        
     end
     
     if ShowFigureFlag
     figure(106);
     subplot(2,2,i*2-1);
     hold on
     %Add scene area
     %Plot the region for scene 
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
        % x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),'-');
         c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
     
     NumEach=cellfun(@numel,EyeChannelS_X_Transform);

     p=cellfun(@plot,EyeChannelS_X_Transform(NumEach>0),EyeChannelS_Y_Transform(NumEach>0));
      for ii=1:length(p)
        % f(ii).Marker='o';
         p(ii).LineStyle='-';
         p(ii).Color=ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:);
         p(ii).Color(4)=0.15;
         p(ii).LineWidth=0.5;
         p(ii).MarkerSize=0.5;
      end
      %Overlay object location
         
         FirstTargetLoc_Aligned_X=0;
         FirstTargetLoc_Aligned_Y=15;
  
         r=7.5;
         t=cplot(r,FirstTargetLoc_Aligned_X,FirstTargetLoc_Aligned_Y);
         t.FaceColor=ColorFractal;
         t.FaceAlpha=0.1;
         t.EdgeColor=ColorFractal;
         t.EdgeAlpha=0.1;
          title(strcat(Blockstr{BlockOrder(i)},',Scramble'),'FontSize',20);
         
         box on
     axis equal;
     axis off
     xlim([-30,30])
     ylim([-30,30])
         
        subplot(2,2,i*2);
     hold on
     %Add scene area
     %Plot the region for scene 
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
        % x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),'-');
         c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
     
     NumEach=cellfun(@numel,EyeChannelN_X_Transform);

     p=cellfun(@plot,EyeChannelN_X_Transform( NumEach>0),EyeChannelN_Y_Transform(NumEach>0));
      for ii=1:length(p)
        % f(ii).Marker='o';
         p(ii).LineStyle='-';
         p(ii).Color=ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:);
         p(ii).Color(4)=0.15;
         p(ii).LineWidth=0.5;
         p(ii).MarkerSize=0.5;
      end
      %Overlay object location
         
         FirstTargetLoc_Aligned_X=0;
         FirstTargetLoc_Aligned_Y=15;
  
         r=7.5;
         t=cplot(r,FirstTargetLoc_Aligned_X,FirstTargetLoc_Aligned_Y);
         t.FaceColor=ColorFractal;
         t.FaceAlpha=0.2;
         t.EdgeColor=ColorFractal;
         t.EdgeAlpha=0.2;  
          title(strcat(Blockstr{BlockOrder(i)},',Normal'),'FontSize',20);
         
    
     
     box on
     axis equal;
     axis off
     xlim([-30,30])
     ylim([-30,30])
     %}
     end
     EyeChannel_X_First_Explore_Block{i}=cellfun(@(x) x{1},EyeChannel_X_Explore,'UniformOutput',0);
     EyeChannel_Y_First_Explore_Block{i}=cellfun(@(x) x{1},EyeChannel_Y_Explore,'UniformOutput',0);
  
end



%Loop over the fractal order for comparison
for i=1:max(TargetNum)
    currTime=cellfun(@(x) x(:,i),SearchingTime_CurrBlock_For_Stats,'uniform',0);
    NumberOfTime=cellfun(@numel,currTime);
    
    currTimeMean=cellfun(@(x) x(:,i),SearchingTime_Mean)';
    currTimeSEM=cellfun(@(x) x(:,i),SearchingTime_SEM)';
    SearchingTime_Max_Order(i)=max(currTimeMean+currTimeSEM);
    
    
    
    
    blocktype_temp=[];
    for j=1:length(NumberOfTime)
        blocktype_temp=[blocktype_temp;UniqueBlockType(j)*ones(NumberOfTime(j),1)];
    end
    currTime_mat=cell2mat(currTime');
    SearchingTime_MaxOrder(i)=max(currTime_mat);
    
    
    
    p_TimeEachFractalOrder(i) = anova1(currTime_mat,blocktype_temp,'off');
    try
    [ p_TimeEachFractalOrder_RankSum(i),~]=ranksum(currTime{ismember(UniqueBlockType,[10,11])},currTime{ismember(UniqueBlockType,[0,1])});
    catch
        p_TimeEachFractalOrder_RankSum(i)=NaN;
    end
    
end


%compare normal with scramble
%STNormal_Mat_Raw=
 for i=1:3
     
     %{
     for j=1:NumberOfBlocks
%for i=3:3
    CurrentBlock=UniqueBlockType(i);
    BlockOrder(i)=find(AssociateBlock==CurrentBlock);
    
     end
     %}
     p_TimeScrambleCompareFractalOrder_Stable(i)=ranksum(STScramble_Mat_Raw{UniqueBlockType==11,i},STNormal_Mat_Raw{UniqueBlockType==11,i});
     p_TimeScrambleCompareFractalOrder_Flexible(i)=ranksum(STScramble_Mat_Raw{UniqueBlockType==1,i},STNormal_Mat_Raw{UniqueBlockType==1,i});
     
 end


%SearchingTime_CurrBlock_For_Stats{i}=


%%%%%%%%%%%%%%%%%%%%%%%%%%%Plottings%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ShowFigureFlag
if ishandle(101) 
    close 101; end
figure(101);set(101,'color','w');
set(101,'Position', [100,100 800,600], 'Name', 'Explore Time as a function of fractal numbers');

%Use barplot
%Reshape Mat
Mat=[];
Mat_tmp=[];
Mat_err=[];
Mat_err_tmp=[]
for i=1:size(STScramble_Mat_Mean,1)
    Mat_tmp=[STScramble_Mat_Mean(i,:)',STNormal_Mat_Mean(i,:)'];
    Mat=[ Mat, Mat_tmp];
    Mat_err_tmp=[STScramble_Mat_SEM(i,:)',STNormal_Mat_SEM(i,:)'];
    Mat_err=[ Mat_err, Mat_err_tmp];
    
end




b=bar(Mat,'grouped');
hold on
index=1;
ngroups = size(Mat, 1);
nbars = size(Mat, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i=1:length(b)
     b(i).FaceColor=ColorForEachMode(BlockOrder(index),:);
    if mod(i,2)==0
  
        index=index+1;

    end
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Mat(:,i), Mat_err(:,i), '.','Color','k','LineWidth',3);
    
end
box off



%}
%{
 STNormal_Mat_Mean(i,iii)=nanmean(SearchingTimeNormal_CurrBlock(:,iii),'all');
             STNormal_Mat_SEM(i,iii)=nanstd(SearchingTimeNormal_CurrBlock(:,iii),'all')/sqrt(sum(~isnan(SearchingTimeNormal_CurrBlock(:,iii)),'all'));
    %         STNormal_Mat_Raw{i,iii}=SearchingTimeNormal_CurrBlock(:,iii);
 %}            

%{
%Add p_value
%p_TimeEachFractalOrder_RankSum(i)
x=1:length(p_TimeEachFractalOrder_RankSum);

yData=SearchingTime_Max_Order*1.05;


text(x(p_TimeEachFractalOrder_RankSum<0.05)-0.045,yData(p_TimeEachFractalOrder_RankSum<0.05)*1.1,'*','FontSize',40,'Color','r');
text(x(p_TimeEachFractalOrder_RankSum<0.01)-0.045,yData(p_TimeEachFractalOrder_RankSum<0.01)*1.2,'*','FontSize',40,'Color','r');
text(x(p_TimeEachFractalOrder_RankSum<0.001)-0.045,yData(p_TimeEachFractalOrder_RankSum<0.001)*1.3,'*','FontSize',40,'Color','r');

ylim([0,max(yData)*1.4]);
%}
%legend(Blockstr( BlockOrder));
xlabel('Fractal Order (N-th)');
ylabel('Mean Exploration TimeSpan(s)');
set(gca,'FontWeight','Bold');
set(gca,'FontSize',25);
set(gca,'LineWidth',3)
xlim([0.5,3+0.5]);



%{
%%
%%%%%%%%%%%%%% Distribution of Exploration TimeSpan 
%%%%%%%%%%%%%% 
titlestr{1}='ExplorationTimeSpanDistribution';
fig{1}=PrepareFigure(107,'w',[120,100, 800,600],'Name',titlestr{1});

EdgesPlot=(Edges(1:end-1)+Bin/2);

for i=1:NumberOfBlocks
     subplot(ceil(NumberOfBlocks/2),2,i);

b=bar(EdgesPlot,SearchingTime_CurrBlock_Count(i,:));
b.FaceColor=ColorForEachMode(BlockOrder(i),:);
b.EdgeColor=b.FaceColor;

%Add mean
hold on
yrange=max(SearchingTime_CurrBlock_Count(i,:),[],'all');
p=plot(SearchingTime_MeanAll(i),yrange*1.1,'V','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');


textstring=sprintf('Mean: %2.1f',SearchingTime_MeanAll(i));
text(SearchingTime_MeanAll(i)*0.8,yrange*1.2,textstring,'FontSize',20,'FontWeight','Bold');
ylim([0,yrange*1.4]);
box off
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
%set(gca,'XScale','log')

if mod(NumberOfBlocks,2)==0
    if i>NumberOfBlocks-2
        xlabel('Exploration Timespan (s)','FontSize',25,'FontWeight','Bold');
       
    end
else
    if i>NumberOfBlocks-1
         xlabel('Exploration Timespan (s)','FontSize',25,'FontWeight','Bold');
       
    end
    
end
if mod(i,2)==1
     ylabel('Number of cases','FontSize',25,'FontWeight','Bold');
    
end
xticks([1:2:SepPoint,SepPoint]);


x_labels={'1','3','5',sprintf('>%d',SepPoint)};
xticklabels(x_labels);


title(Blockstr{BlockOrder(i)},'FontSize',25,'FontWeight','Bold');
    
end

%%

%%%%%%%%%%%%%% Overlap Distribution of Exploration TimeSpan 
%%%%%%%%%%%%%% 
titlestr{2}='ExplorationTimeSpanDistributionCompare';
fig{2}=PrepareFigure(108,'w',[140,100, 800,600],'Name',titlestr{2});

b=bar(EdgesPlot',(SearchingTime_CurrBlock_Count./sum(SearchingTime_CurrBlock_Count,2))','grouped');
for i=1:NumberOfBlocks
    b(i).FaceColor=ColorForEachMode(BlockOrder(i),:);
    b(i).EdgeColor=b(i).FaceColor;
    b(i).BarWidth=0.8;
    %Add mean
    hold on
    yrange=max(SearchingTime_CurrBlock_Count./sum(SearchingTime_CurrBlock_Count,2),[],'all');

    plot(SearchingTime_MeanAll(i)',yrange,'V','MarkerSize',10,'MarkerFaceColor',ColorForEachMode(BlockOrder(i),:),'MarkerEdgeColor',ColorForEachMode(BlockOrder(i),:));
     textstring=sprintf('Mean:\n %2.1f',SearchingTime_MeanAll(i));
     text(SearchingTime_MeanAll(i)-0.2,yrange*1.2,textstring,'FontSize',25,'FontWeight','Bold','Color',ColorForEachMode(BlockOrder(i),:));
     
end
ylim([0,yrange*1.4]);
yticks([0:0.5:yrange*1.4]);
     box off
     set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
     xlabel('Exploration Timespan (s)','FontSize',25,'FontWeight','Bold');
     ylabel('Proportion of cases','FontSize',25,'FontWeight','Bold');
     legend(Blockstr{BlockOrder});
     
%b.FaceColor=ColorForEachMode(BlockOrder(i),:);
%b.EdgeColor=b.FaceColor;



%%

%%%%%%%%%%%%%% Overlap Distribution of Exploration TimeSpan 
%%%%%%%%%%%%%% 
if length(SearchingTime_CurrBlock_Line)>1
titlestr{3}='ExplorationTimeSpanMeanSEM';
fig{3}=PrepareFigure(109,'w',[150,100, 800,600],'Name',titlestr{3});

%Temporally put here
%[h,p]=ttest2(SearchingTime_CurrBlock_Line{1}',SearchingTime_CurrBlock_Line{2}');
[h,p]=ranksum(SearchingTime_CurrBlock_Line{1}',SearchingTime_CurrBlock_Line{2}');

b=bar([1:length(SearchingTime_MeanAll)],SearchingTime_MeanAll);
b.FaceColor='Flat';
b.CData=ColorForEachMode(BlockOrder,:);
b.EdgeColor='Flat';


hold on
errorbar([1:length(SearchingTime_MeanAll)],SearchingTime_MeanAll,SearchingTime_SemAll,'sk','LineWidth',3);
 box off
 set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
    
ylabel('Mean & SEM','FontSize',25,'FontWeight','Bold');
 xticklabels([Blockstr{BlockOrder}]);
 yrange=ylim;
 %Add significance level
 if p<0.05
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on

yrange=ylim;
plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',3)
if p<0.05 & p>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.2, '*k','MarkerSize',20,'LineWidth',3)
elseif p<0.01 & p>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.2, '*r','MarkerSize',20,'LineWidth',3)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.2, '*r','MarkerSize',20,'LineWidth',3)
        
elseif p<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.2, '*r','MarkerSize',20,'LineWidth',3)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.2, '*r','MarkerSize',20,'LineWidth',3)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.2, '*r','MarkerSize',20,'LineWidth',3)
    
end
 end
set(gca,'xticklabel',{[]});
end


%%%%%%%%%%%%%% Distribution of TargetNum
%%%%%%%%%%%%%% 
titlestr{4}='Distribution of target num';
fig{4}=PrepareFigure(110,'w',[170,100, 800,600],'Name',titlestr{4});
h=histogram(TargetNum,'FaceColor','k','EdgeColor','k');
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
xlabel('Number of fractals appeared in each trial','FontSize',30);
ylabel('Number of cases','FontSize',30);
box off
%}
end %End of the if ShowFigureFlag


     


end