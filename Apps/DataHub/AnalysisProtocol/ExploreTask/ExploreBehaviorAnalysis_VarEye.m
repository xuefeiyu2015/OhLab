function  ExploreBehaviorAnalysis_VarEye(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Xuefei Yu, 03,21,2020
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

%FP Location
FixLoc_X=squeeze(Data.FixLoc(Selection,2:end,1));
FixLoc_Y=squeeze(Data.FixLoc(Selection,2:end,2));

[FPAngle_tmp,FPEcc_tmp]=LocToAngle(FixLoc_X,FixLoc_Y);
%In this task, the FP angle and ecc are fix for each trial
%{
FPAngle=FPAngle_tmp(:,1);
FPEcc=FPEcc_tmp(:,1);
%}
FPAngle=FPAngle_tmp;
FPEcc=FPEcc_tmp;


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
 %{
 FPAngle=FPAngle(AdvanceInTargetTrial);
FPEcc=FPEcc(AdvanceInTargetTrial);
 %}
 FPAngle=FPAngle(AdvanceInTargetTrial,:);
 FPEcc=FPEcc(AdvanceInTargetTrial,:);
 
FixLoc_X=FixLoc_X(AdvanceInTargetTrial,:);
FixLoc_Y=FixLoc_Y(AdvanceInTargetTrial,:);

Tgt_X=TargetLoc(:,:,1);
Tgt_Y=TargetLoc(:,:,2);

DistFPToTarget=sqrt((Tgt_X-FixLoc_X).^2+(Tgt_Y-FixLoc_Y).^2);
DistFPToTargetRef=uniquetol(DistFPToTarget(FPEcc==0),0.001);
RelDistFPToTarget=DistFPToTarget./DistFPToTargetRef;

 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime)/1000;%Change into s for presentation

RelSearchingTime=SearchingTime./RelDistFPToTarget;


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

%Loop over the FP eccentricity and blocks
FPEcc_First=FPEcc(:,1);
FPAngle_First=FPAngle(:,1);
UniqueFPEcc=unique(FPEcc_First);
UniqueAngle=unique(FPAngle_First(FPEcc_First~=0));
for i=1:length(UniqueFPEcc)
    
    for j=1:NumberOfBlocks
        
            Sel=FPEcc(:,1)==UniqueFPEcc(i) & BlockType==UniqueBlockType(j);
            SearchingTime_Sel=RelSearchingTime(Sel,:);
            %Loop over the fractal order and statistics
            %Only use order1,2,>2
            for m=1:3
                if m<3
                  SearchingTime_Sel_Mean{i,j}(m)=nanmean(SearchingTime_Sel(:,m));
                   SearchingTime_Sel_Sem{i,j}(m)=nanstd(SearchingTime_Sel(:,m))./sqrt(sum(~isnan(SearchingTime_Sel(:,m))));
                   SearchngTime_Ecc{i,j,m}=RelSearchingTime(Sel,m);
                else
                  SearchingTime_Sel_Mean{i,j}(m)=nanmean(SearchingTime_Sel(:,m:end),'all');
                  SearchingTime_Sel_Sem{i,j}(m)=nanstd(SearchingTime_Sel(:,m:end),[],'all')./sqrt(sum(~isnan(SearchingTime_Sel(:,m:end)),'all'));
                  SearchngTime_Ecc{i,j,m}=reshape(RelSearchingTime(Sel,m:end),numel(RelSearchingTime(Sel,m:end)),1); 
                end
               
              

            end

            
            
            
            
           
        %{
        if UniqueFPEcc(i)~=0
         for m=1:length(UniqueAngle)
        
            Sel=FPEcc(:,1)==UniqueFPEcc(i) & BlockType==UniqueBlockType(j) & FPAngle(:,1)==UniqueAngle(m);
            SearchingTime_Sel=RelSearchingTime(Sel,:);
            SearchingTime_Sel_PMean{i,j,m}=nanmean(SearchingTime_Sel,1);
            SearchingTime_Sel_PSem{i,j,m}=nanstd(SearchingTime_Sel,[],1)./sqrt(sum(~isnan(SearchingTime_Sel),1));
        end
        end
 %}
    end

end

%Loop over the fractal order for statistics
%Only use order1,2,>2

for j=1:length(UniqueFPEcc)
 for i=1:3
    try
    [p_FractalOrder(j,i),h]=ranksum(SearchngTime_Ecc{j,1,i},SearchngTime_Ecc{j,2,i});
    catch
    p_FractalOrder(j,i)=NaN; 
    end
 end
end


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

%Next Eye trace align with each target

for i=1:NumberOfBlocks
%for i=3:3
    CurrentBlock=UniqueBlockType(i);
    BlockOrder(i)=find(AssociateBlock==CurrentBlock);
   
    
     %Select out the eye movements during the searching period for each
     %scene 
     SceneName_CurrBlock=SceneName(BlockType==CurrentBlock,:);
     UniqueScene_CurrBlock=unique(SceneName_CurrBlock);
     EyeChannel_X_CurrBlock=EyeChannel_X_L(BlockType==CurrentBlock);
     EyeChannel_Y_CurrBlock=EyeChannel_Y_L(BlockType==CurrentBlock);
     
     %Select out the object location
     TargetLoc_CurrBlock=TargetLoc(BlockType==CurrentBlock,:,:);
     TargetPolarAngle_CurrBlock=TargetPolarAngle(BlockType==CurrentBlock,:,:);
     
     %Select out the fp location
     FP_XCurrBlock=FixLoc_X(BlockType==CurrentBlock,:);
     FP_YCurrBlock=FixLoc_Y(BlockType==CurrentBlock,:);
     
     FP_EccCurrBlock=FPEcc(BlockType==CurrentBlock,:,:);
     FP_AngleCurrBlock=FPAngle(BlockType==CurrentBlock,:,:);
  
     FixOffTime_CurrBlock=FixOffTime(BlockType==CurrentBlock,:);
     TargetOnTime_CurrBlock=TargetOnTime(BlockType==CurrentBlock,:);
     
     
     
     for j=1:length(UniqueScene_CurrBlock)
     % for j=3:3
         %Select the eye trace for each scene
         %Remove ecc==0
         Select_Scene=(SceneName_CurrBlock==UniqueScene_CurrBlock(j))&(FP_EccCurrBlock(:,1)~=0);
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
         %Select out current FP location
         FirstFPLoc_CurrScene=[FP_XCurrBlock(Select_Scene,FractalOrder),FP_YCurrBlock(Select_Scene,FractalOrder)];
         FPPolarAngleEcc_CurrScene=[FPAngle(Select_Scene,FractalOrder),FPEcc(Select_Scene,FractalOrder)];
         if size(FirstFPLoc_CurrScene,2)==1
             FirstFPLoc_CurrScene=FirstFPLoc_CurrScene';
         end
         
         %Plot, put it here to avoid another loop
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
        
        %
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
         %}
         %Overlay object location
       
         FirstTargetLoc_CurrSceneUnique=uniquetol(FirstTargetLoc_CurrScene,'ByRows',true);
         
         FirstTargetLoc_CurrScene_X=FirstTargetLoc_CurrSceneUnique(:,1);
         FirstTargetLoc_CurrScene_Y=FirstTargetLoc_CurrSceneUnique(:,2);
         
         
         
         for iii=1:length(FirstTargetLoc_CurrScene_X)
             if ~isnan(FirstTargetLoc_CurrScene_X(iii))
             
                r=5;
                t=cplot(r,FirstTargetLoc_CurrScene_X(iii),FirstTargetLoc_CurrScene_Y(iii));
                t.FaceColor=ColorFractal;
                t.FaceAlpha=0.2;
                t.EdgeColor=ColorFractal;
                t.EdgeAlpha=0.1 ;
             end
             
         end
         %Overlay fixation point
         FirstFPLoc_CurrSceneUnique=uniquetol(FirstFPLoc_CurrScene,'ByRows',true);
         FirstFPLoc_CurrScene_X=FirstFPLoc_CurrSceneUnique(:,1);
         FirstFPLoc_CurrScene_Y=FirstFPLoc_CurrSceneUnique(:,2);
          for iii=1:length(FirstFPLoc_CurrSceneUnique)
             if ~isnan(FirstFPLoc_CurrSceneUnique(iii))
             
                r=1;
                t=cplot(r,FirstFPLoc_CurrScene_X(iii),FirstFPLoc_CurrScene_Y(iii));
                t.FaceColor='#FB9F41';
               % t.FaceAlpha=0.1;
                t.EdgeColor='#FB9F41';
                %t.EdgeAlpha=1;
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
       temp_AfterTransform_X=[];
       temp_AfterTransform_Y=[];
       
      
       
       for iii=1:length(TargetPolarAngle_CurrScene)
          if ~isnan(TargetPolarAngle_CurrScene(iii))
           TransformMat=[cos(TargetPolarAngle_CurrScene(iii)*pi/180),-sin(TargetPolarAngle_CurrScene(iii)*pi/180);sin(TargetPolarAngle_CurrScene(iii)*pi/180),cos(TargetPolarAngle_CurrScene(iii)*pi/180)];
       
           EyeChannel_X_First_Explore_Before=EyeChannel_X_First_Explore{iii};
           EyeChannel_Y_First_Explore_Before=EyeChannel_Y_First_Explore{iii};
           if isempty(EyeChannel_X_First_Explore_Before)
               temp_AfterTransform_X=[NaN];
               temp_AfterTransform_Y=[NaN];
           else
               %{
                temp=[(EyeChannel_X_First_Explore_Before-FirstFPLoc_CurrScene(iii,1))',(EyeChannel_Y_First_Explore_Before-FirstFPLoc_CurrScene(iii,2))']*TransformMat;
                 temp_AfterTransform_X=[temp_AfterTransform_X+FirstFPLoc_CurrScene(iii,1),temp(:,1)',NaN];
                temp_AfterTransform_Y=[temp_AfterTransform_Y+FirstFPLoc_CurrScene(iii,2),temp(:,2)',NaN];
               %}
                temp=[(EyeChannel_X_First_Explore_Before)',(EyeChannel_Y_First_Explore_Before)']*TransformMat;
                 temp_AfterTransform_X=[temp_AfterTransform_X,temp(:,1)',NaN];
                temp_AfterTransform_Y=[temp_AfterTransform_Y,temp(:,2)',NaN];
               
           
           end
       
    
          end
       
       end
       
       
       EyeChannel_X_Transform{j}=temp_AfterTransform_X;
       EyeChannel_Y_Transform{j}=temp_AfterTransform_Y;
         
         
     
         
         
     end
     
      if ShowFigureFlag
     figure(106);
     subplot(2,2,i);
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
     
     

     p=cellfun(@plot,EyeChannel_X_Transform,EyeChannel_Y_Transform);
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
  
         r=5;
         t=cplot(r,FirstTargetLoc_Aligned_X,FirstTargetLoc_Aligned_Y);
         t.FaceColor=ColorFractal;
         t.FaceAlpha=0.2;
         t.EdgeColor=ColorFractal;
         t.EdgeAlpha=0.2;
             
         
     
      
      title(Blockstr{BlockOrder(i)},'FontSize',20);
     
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






%%Figure
FigureStartNum=110;
FigureIndex=1;

if ShowFigureFlag
 %Figure 1: 
figtitlestr{FigureIndex}='RTWithEccentricity';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});    

for i=1:length(UniqueFPEcc)
    subplot(length(UniqueFPEcc),1,i)
    
    for j=1:NumberOfBlocks
        
       l=errorbar([1:length(SearchingTime_Sel_Mean{i,j})],SearchingTime_Sel_Mean{i,j}, SearchingTime_Sel_Sem{i,j},'Color',ColorForEachMode(find(AssociateBlock==UniqueBlockType(j)),:),'Marker','s','MarkerSize',15);
       
    
    l.MarkerEdgeColor=ColorForEachMode(find(AssociateBlock==UniqueBlockType(j)),:);
    l.MarkerFaceColor=l.MarkerEdgeColor;
   
    l.LineWidth=3;
       hold on
       %p_FractalOrder(j,i)
       title(sprintf('FP Eccentricity %d',UniqueFPEcc(i)));
       ylabel('Relative RT(a.u.)');
       if i==length(UniqueFPEcc)
           xlabel('Fractal Order(N-th)');
       end
       set(gca,'FontWeight','Bold');  
set(gca,'FontSize',15);
set(gca,'LineWidth',3);
xlim([0.5,3+0.5]);
box off;
 
       
    end
    %Statistics
    x=1:length(p_FractalOrder(i,:));
     maxy=ylim;
    yData=maxy(2)*1.05*ones(1,length(x));


    text(x(p_FractalOrder(i,:)<0.05),yData(p_FractalOrder(i,:)<0.05)*1.1,'*','FontSize',40,'Color','r');
    text(x(p_FractalOrder(i,:)<0.01),yData(p_FractalOrder(i,:)<0.01)*1.2,'*','FontSize',40,'Color','r');
    text(x(p_FractalOrder(i,:)<0.001),yData(p_FractalOrder(i,:)<0.001)*1.3,'*','FontSize',40,'Color','r');

    ylim([0,max(yData)*1.4]);


end

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;
%{
%Figure 2: 
figtitlestr{FigureIndex}='RTWithEachAngle';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 2000,800],'Name',figtitlestr{FigureIndex});    

for i=1:length(UniqueAngle)
    subplot(length(UniqueAngle),1,i)
    
    for j=1:NumberOfBlocks
        
       errorbar([1:length(SearchingTime_Sel_PMean{2,j,i})],SearchingTime_Sel_PMean{2,j,i}, SearchingTime_Sel_PSem{2,j,i},'Color',ColorForEachMode(find(AssociateBlock==UniqueBlockType(j)),:),'Marker','s','MarkerSize',15);
       hold on
    end

end


FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;
%}
end %End of the if figure plot flag

end


