function VideoReconstructionForPath(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

%This is to reconstruct monkey's performance during path tasks
%task;
%Xuefei Yu; 01182025
%Unpack all the data base and task settings
ECode;


%{
%Panel control mode
%
if EndTrial-StartTrial<=8
    ReConstructTrial=[StartTrial: EndTrial];
else
    ReConstructTrial=[EndTrial-4: EndTrial];
end
%}
%ReConstructTrial=[13,29,37,15,83,84,97,101];
ReConstructTrial=[StartTrial,EndTrial];

NeutralFP=1;%Use Neutral FP

%Before learning v.s. After Learning Trials
 LearningPlot=1;
% LearningStr='Late Stage of Training......';
 LearningStr='After 30 days......';
% LearningStr='Initial Learning';

%ReConstructTrial=[906,382,[FixLocTrials(5:6)]];%[906];%[FixLocTrials(5:6)];
%FixLocTrials=[351,339,353,20];%[351,353,339];% For the demonstration of
%performance before and after training.s
%ReConstructTrial=FixLocTrials;

%{
%Create a file to save video
writerObj=VideoWriter('test.mp4','MPEG-4');  
SaveVideo=1;
if SaveVideo
    open(writerObj);
end
%}

TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel
SaveFlag=0;%Flag to indicate whether saving the data

%Hard wired parameters
TargetOnOffset=50;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

ScreenWidth=80; %in degree
ScreenHeight=80;%in degree

SceneAreaWidth=50;%in degree
SceneAreaHeight=50;%in degree

FractalWidth=10; %in degree
FractalHeight=10; %in degree

FixationWindow=8; %in degree
TargetWindow=15; %Diameter; in degree

RewardAmount=200; %in ms




%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');
%To get the monkey name
RemainFiles=erase(Data.FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);





%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(:,1);
EyeChannel_Y_L=EyeData(:,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

EyeBinWidth=Data.EyeBinWidth;

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);
%Ensure that the eventbin start from 1
DiffFirst=EventBin(:,1)-1;
EventBin=EventBin-DiffFirst;




EyeChannel_X_L=EyeChannel_X_L(Selection);
EyeChannel_Y_L=EyeChannel_Y_L(Selection);

%Get scene and fractals information
%SceneName=Data.SceneIndex(Selection,3);
SceneName=1;
TargetName=Data.TargetIndex(Selection,2:end,2);
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%Sequence information from the library
ParaLib=Data.ParaLib;
if ismember('HyperSetNO,ManuelFlag,LocSeqCode',keys(ParaLib))
    Seq_Info=ParaLib('HyperSetNO,ManuelFlag,LocSeqCode');
    Seq_Info=Seq_Info(Selection,:);
    SeqNo=Seq_Info(:,2);
end
SeqNoDisplay = unique(SeqNo)

ExampleSeqNo = 19;

CandidateTrials = find(SeqNo==ExampleSeqNo);


TargetNum=Data.TargetNum(Selection,2:end);
TargetSize=Data.TargetSize(Selection,2:end);
TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);
TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);

%Set up the path for scene and fractals

Workingdirectory=pwd;

%Set up output path according to the system
OperationSystem = computer;%Get the operation system information

if strcmp(OperationSystem(1:3),"PCW")  
    %PC
     MarkerFolder='DataHub';
     Flag=strfind(Workingdirectory,MarkerFolder);
     BasicDirectory=Workingdirectory(1:Flag-1);
elseif strcmp(OperationSystem(1:3),"MAC")  
    %Mac
     MarkerFolder='DataAnalysis';
     Flag=strfind(Workingdirectory,MarkerFolder);
     BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));
 end



ImageFolerName='VisualStimulus';
ImagePath=strcat(BasicDirectory,ImageFolerName,BasicDirectory(end));
%{
%Select trials for reconstruction
try
TrialsForExplore=find(BlockType==1|BlockType==0);
TrialsForFixLoc=find(BlockType==11|BlockType==10);

FractalNumForExplore=TargetNum(TrialsForExplore);
FractalNumForFixLoc=TargetNum(TrialsForFixLoc);


AdvanceInTargetNum=sum(EventChannel==ADVEYE_IN,2);
AdvanceInTargetTrial=find(AdvanceInTargetNum==TargetNum);

%Select Searching time
 FixOffTime=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);
 TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
 
SearchingTime= TargetOnTime-FixOffTime;
SearchingTimeForExplore=SearchingTime(TrialsForExplore,:);

%SortTrialAccordingToSearchingTime

[SortRearchingTime,Index]=sort(SearchingTimeForExplore(:,1));
Trial_Explore_Sort=TrialsForExplore(Index);

MiddleTrial=Trial_Explore_Sort(ceil(length(Trial_Explore_Sort)/2));
TrialsForFractalNumMoreThanOne=TrialsForExplore(FractalNumForExplore>1);

%Information for selected trials for reconstruction
BlockTypeForRandLoc=BlockType(RandLocTrials);
SceneNameForRandLoc_Selection=SceneName(RandLocTrials);



%Find out the trials with the same scene
TrialsFixLocWithSameScene=find((ismember(SceneName,SceneNameForRandLoc_Selection)&(BlockType==11|BlockType==10))==1);

SceneNameForFixLoc_Selection=SceneName(TrialsFixLocWithSameScene);
FixLocTrialsForFractalNumMoreThanOne=TrialsForFixLoc(FractalNumForFixLoc>1);

TrialsFixLocSelection_MoreThanOne=FixLocTrialsForFractalNumMoreThanOne(ismember(FixLocTrialsForFractalNumMoreThanOne,TrialsFixLocWithSameScene));
catch
end
%SearchingTime_ForCandidateTrial=SearchingTime(TrialsFixLocWithSameScene,1);

SearchingTime_ForCandidateTrial=SearchingTime(TrialsForFixLoc,1);
SearchingTime_ForCandidateTrialAll=SearchingTime(TrialsForFixLoc,:);

[SortRearchingTimeFix,Index]=sort(SearchingTime_ForCandidateTrial,'descend' );
%Trial_Fix_Sort=TrialsFixLocWithSameScene(Index);
Trial_Fix_Sort=TrialsForFixLoc(Index);

%{
%Manually find one trials
SceneNameForSortedTrials=SceneName(Trial_Fix_Sort);
SelectTrial=Trial_Fix_Sort(SceneNameForSortedTrials==1057);
SelectTrial_Advance=SelectTrial(ismember(SelectTrial,AdvanceInTargetTrial));

Trial_Fix_SelectionInfo=SceneName(FixLocTrials);

%}

%}


%Set up the screen patch
if ishandle(100) 
    close(100); end
figure(100);set(100,'color','k');
set(100,'Position', [100,100 600,600], 'Name', 'Path Sequence Task');


%Set up the screen
ax=axes('pos',[0,0,1,1]);
ax.Color=[0,0,0];
axis([-ScreenWidth/2,ScreenHeight/2,-ScreenWidth/2,ScreenHeight/2]);

axis equal;



    


if  LearningPlot
       
        pause(2);
        Learning_Block= text(0,0,LearningStr,'FontWeight','bold','FontSize',25,'Color','w','HorizontalAlignment','center');
        pause(2);
       Learning_Block.Visible='off';
end


%Main reconstruction procedural;
for i=1:length(ReConstructTrial)
   
   % RTrial=ReConstructTrial(i);
   RTrial=i;
   % BlockType_curr=BlockType(RTrial);
    
       
  
    %{
    if (BlockType_curr==0||BlockType_curr==1)
        BlockStr1='Treasure Seeking Task';
        BlockStr2='Flexible Location Series Scene Set Block';
        if ~ShowRandLocBlockStr
        t_Block1= text(0,10,BlockStr1,'FontWeight','bold','FontSize',25,'Color','w','HorizontalAlignment','center');
        t_Block2= text(0,-2,BlockStr2,'FontWeight','bold','FontSize',20,'Color','w','HorizontalAlignment','center');
        
        pause(2);
         t_Block1.Visible='off';
         t_Block2.Visible='off';
         
         pause(1);
        
        
        ShowRandLocBlockStr=1;
        
        end
    else
        BlockStr1='Treasure Seeking Task';
        BlockStr2='Stable Location Series Scene Set Block';
        if ~ShowFixBlockStr
       pause(2);
        t_Block1= text(0,10,BlockStr1,'FontWeight','bold','FontSize',25,'Color','w','HorizontalAlignment','center');
        t_Block2= text(0,-2,BlockStr2,'FontWeight','bold','FontSize',20,'Color','w','HorizontalAlignment','center');
        ShowFixBlockStr=1;
         pause(2);
         t_Block1.Visible='off';
         t_Block2.Visible='off';
         
         pause(1);
        
        
        end
 
    end
    
    %}
    %{
    %Trial Information
     if abs(BlockType_last-BlockType_curr)>1
         TrialIndicator=1;
     else
         TrialIndicator=TrialIndicator+1;
          
         
     end
    %}
   TrialIndexStr=sprintf('#Trial %d',RTrial);
    t= text(0,0.5,TrialIndexStr,'FontSize',40,'Color','w','HorizontalAlignment','center');
    pause(1);
    t.Visible='off';
    pause(0.5);
    
    
    
    
   %
   SceneName_CurrTrial=SceneName;%(RTrial,:);
  % SceneName_CurrTrial=SceneName(i,:);
    SceneNameStr=sprintf('Scene%d.gif',SceneName_CurrTrial);
    scene=imread(strcat(ImagePath,SceneNameStr));
    SceneHight=size(scene,1);
    SceneWidth=size(scene,2);
   %}
    %
    %PreLoad scene patch
    ScenePatch=axes('pos',[(ScreenWidth-SceneAreaWidth)/ScreenWidth/2,(ScreenWidth-SceneAreaWidth)/ScreenHeight/2,SceneAreaWidth/ScreenWidth,SceneAreaWidth/ScreenHeight]);
    ScenePatch.Color='None';
    axis off;  
    SceneImage=imshow(scene);
    SceneImage.Visible='off';
   
    hold on
    %}
    %Preload Fixation point, but make it hide 
   % if NeutralFP
         %Use neutral fixation point
         FPimage=plot(SceneHight/2,SceneWidth/2,'o','MarkerFaceColor','#FBE041','MarkerEdgeColor','#FBE041','MarkerSize',20,'Visible','off');
    %{    
    else
         %Use distint FP Point
        if BlockType_curr==10 || BlockType_curr==11
        FPimage=plot(SceneHight/2,SceneWidth/2,'s','MarkerFaceColor','#FFB000','MarkerEdgeColor','#FFB000','MarkerSize',20,'Visible','off');
        else
        FPimage=plot(SceneHight/2,SceneWidth/2,'+','MarkerFaceColor','#1AFF1A','MarkerEdgeColor','#1AFF1A','MarkerSize',20,'LineWidth',50,'Visible','off');
        end
    end
    %}
   
    
   %
    %Preload Extra Fixation point
    %Get the location of the extra fixation point;
    ExtraFPLoc_CurrTrial=[0,0];    
    ExtraFPLoc_CurrTrial_XRel=ExtraFPLoc_CurrTrial(1)/SceneAreaWidth*SceneWidth;
    ExtraFPLoc_CurrTrial_YRel=-ExtraFPLoc_CurrTrial(2)/SceneAreaWidth*SceneWidth;
     
    ExtraFPimage=plot(SceneHight/2+ExtraFPLoc_CurrTrial_XRel,SceneWidth/2+ExtraFPLoc_CurrTrial_YRel,'o','MarkerEdgeColor','#1A85FF','MarkerSize',20,'LineWidth',5,'Visible','off');
   
   %}
     %Preload fixation window
      FixationWidowRel=FixationWindow/SceneAreaWidth*SceneWidth;
      FixWindow=rectangle('Position',[SceneHight/2-FixationWidowRel/2 SceneHight/2-FixationWidowRel/2 FixationWidowRel FixationWidowRel],...
       'EdgeColor','#FF1A1A','LineWidth',2,'Visible','off');
   %
     %Preload Extra fixation window
     
     
     ExtraFixWindow=rectangle('Position',[ExtraFPLoc_CurrTrial_XRel+SceneHight/2-FixationWidowRel/2 ExtraFPLoc_CurrTrial_YRel+SceneHight/2-FixationWidowRel/2 FixationWidowRel FixationWidowRel],...
         'EdgeColor','#FF1A1A','LineWidth',2,'Visible','off');
   %}
     
    
    %Preload all Fractals
    sel_target = ~isnan(TargetNum(i,:));
    TargetNum_CurrTrial=sum(sel_target);
    TargetName_CurrTrial=TargetName(i,1:TargetNum_CurrTrial*2);
    
    TargetLoc_CurrTrial=squeeze(TargetLoc(i,1:TargetNum_CurrTrial*2,:));
    
     %Preload Target window
       for j=1:TargetNum_CurrTrial*2
         TargetWidowRel=TargetWindow/SceneAreaWidth*SceneWidth;
         TargetLoc_CurrTrial_XRel=TargetLoc_CurrTrial(j,1)/SceneAreaWidth*SceneWidth;
         TargetLoc_CurrTrial_YRel=-TargetLoc_CurrTrial(j,2)/SceneAreaWidth*SceneWidth;

         
         if mod(j,2)==0
             %Error choice
              TgtWindow{j}=rectangle('Position',[TargetLoc_CurrTrial_XRel+SceneHight/2-TargetWidowRel/2 TargetLoc_CurrTrial_YRel+SceneHight/2-TargetWidowRel/2 TargetWidowRel TargetWidowRel],...
         'Curvature',[1 1],'EdgeColor','#00CC66','LineWidth',2,'Visible','off');
        
         else
             %Correct choice
         TgtWindow{j}=rectangle('Position',[TargetLoc_CurrTrial_XRel+SceneHight/2-TargetWidowRel/2 TargetLoc_CurrTrial_YRel+SceneHight/2-TargetWidowRel/2 TargetWidowRel TargetWidowRel],...
         'Curvature',[1 1],'EdgeColor','#FFB000','LineWidth',2,'Visible','off');
         end
      %Preload target cue 
    % TgtCue{j}=plot(TargetLoc_CurrTrial_XRel+SceneHight/2,TargetLoc_CurrTrial_YRel+SceneHight/2 ,'s','MarkerFaceColor','#41A2FB','MarkerEdgeColor','#41A2FB','MarkerSize',20,'Visible','off');
       
       
       end
    
    
    
    
    
    for j=1:TargetNum_CurrTrial*2
        
         TargetName_Order=TargetName_CurrTrial(j);
         TargetLoc_Order=TargetLoc_CurrTrial(j,:);
         %Retrive real name
         if TargetName_Order<10000
            SetPart=floor(TargetName_Order/10);
           NoPart=mod(TargetName_Order,SetPart);
 
         else
           SetPart=floor(TargetName_Order/100);
           NoPart=mod(TargetName_Order,SetPart);
        
        
         end
   
         
         TargetNameStr=sprintf('Fractal%d_%d.png',SetPart,NoPart);
         [Fractal,~,alpha]=imread(strcat(ImagePath,TargetNameStr));
          fractal_ax=axes('pos',[(ScreenWidth-FractalWidth)/ScreenWidth/2+TargetLoc_Order(1)/ScreenWidth,(ScreenWidth-FractalWidth)/ScreenHeight/2+TargetLoc_Order(2)/ScreenWidth,...
              FractalWidth/ScreenWidth,FractalWidth/ScreenHeight]);
           fractal_ax.Color='None';
          axis off;
         f{j} = imshow(Fractal);
         set(f{j}, 'AlphaData', alpha,'Visible','off');
   
        
    end
    
  
    
    
    
    %Load event and event time for current trial
    Event_CurrTrial=EventChannel(i,:);%Get all events for current trial
    EventBinChannel_CurrTrial=EventBin(i,:);%Get all events time for current trial
       
    EyeChannel_X_CurrTrial=EyeChannel_X_L{i};
    EyeChannel_Y_CurrTrial=EyeChannel_Y_L{i}+2;
    
 
    EraseIndex=1;
    EraseLength=500;
    ExtraFPIndex = 0;
    %Set up the screen
    axplot=axes('pos',[0,0,1,1]);
    axplot.Color='None';
    axis([-40,40,-40,40]);
    h = animatedline('Marker','o','Color','r');
    h_eye = animatedline('Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',15);
 
    EventPointer=1;
    TargetCounter=1;
    for ii=1:length(EyeChannel_X_CurrTrial)
        %Plot lively eye trace
        hold on
        %Check major event
        if ismember(ii,EventBinChannel_CurrTrial(EventPointer:end))
            
            Event_CurrSelect=Event_CurrTrial(EventPointer:end);
            Selector=find(EventBinChannel_CurrTrial(EventPointer:end)==ii);
            for j=1:length(Selector)
                
                 CurrentEvent=Event_CurrSelect(Selector(j));
     
           
               
            if CurrentEvent==SHOWFIXCD
                           
               FPimage.Visible='on';
               FixWindow.Visible='on'; 
               % clearpoints(h);
               % EraseIndex=0;
               %{
            elseif CurrentEvent==PREFIX
                 FixWindow.Visible='on';            
                 FPimage.Visible='on';
                
               %}
            elseif  CurrentEvent==FIXOFF
                FixWindow.Visible='off';
                 FPimage.Visible='off';

                 if ExtraFPIndex==1
                     ExtraFixWindow.Visible='off';
                 ExtraFPimage.Visible='off';
                 end
                 
                 clearpoints(h);
                EraseIndex=0;
                %{
                TgtWindow{TargetCounter}.Visible='on';
                TgtWindow{TargetCounter+1}.Visible='on';
                %}
                
            elseif CurrentEvent==SCENE_ON  
                SceneImage.Visible='on';
                 clearpoints(h);
                EraseIndex=0;
                
                
            elseif CurrentEvent==SCENE_OFF
                SceneImage.Visible='off';
                 clearpoints(h);
                EraseIndex=0;
                %{
            elseif CurrentEvent==CUEON  
                TgtCue{TargetCounter}.Visible='on';
                clearpoints(h);
                EraseIndex=0;
                %}
            elseif CurrentEvent>=FIRSTTGT & CurrentEvent<=FIFTHTGT
                 f{TargetCounter}.Visible='on';
                 f{TargetCounter+1}.Visible='on';
                 TgtWindow{TargetCounter}.Visible='on';
              %   TgtWindow{TargetCounter+1}.Visible='on';

               

                
                 % clearpoints(h);
               % EraseIndex=0;

               %
           % elseif CurrentEvent==CHOICECONFIRM
%
            elseif CurrentEvent==TGTOFF 
                 TargetCounter=TargetCounter+2;
                %{
                if TargetCounter~=TargetNum_CurrTrial
                    %FP window on
                    FixWindow.Visible='on';
                end
                %}
                %{
                f{TargetCounter}.Visible='off';
                 TgtWindow{TargetCounter}.Visible='off';

                 f{TargetCounter+1}.Visible='off';
                 TgtWindow{TargetCounter+1}.Visible='off';
                
                TargetCounter=TargetCounter+2;
                 clearpoints(h);
%}
              %  EraseIndex=0;
               %}
              %{
            elseif CurrentEvent==REJECT || CurrentEvent==ACCEPT
                TgtWindow{TargetCounter}.Visible='off';
              %}
            elseif CurrentEvent == FINALFIX
                 ExtraFPimage.Visible='on';
                 ExtraFixWindow.Visible='on'; 

                 ExtraFPIndex = 1;

            elseif CurrentEvent==RWDCD
                if ExtraFPIndex==1
                    %pause(0.1);
                    Beeper(800,1,0.2);
                   % pause(0.2);
                else
                   % pause(0.1);
                    Beeper(800,1,0.04);
                  %  pause(0.2);
                end
                 f{TargetCounter}.Visible='off';
                 TgtWindow{TargetCounter}.Visible='off';

                 f{TargetCounter+1}.Visible='off';
                % TgtWindow{TargetCounter+1}.Visible='off';
                
               
                 clearpoints(h);
               % Beeper(1000,100,10);
              
            elseif CurrentEvent==ERR_BEEP
               % pause(0.1);
                Beeper(300,1,0.2);
                 f{TargetCounter}.Visible='off';
                 TgtWindow{TargetCounter}.Visible='off';

                 f{TargetCounter+1}.Visible='off';
             %   pause(0.2);


                
            end
            
            end
            
            
            
            
            
            
            
            
            EventPointer=EventPointer+1;
            
        end
        
  
        addpoints(h,EyeChannel_X_CurrTrial(ii),EyeChannel_Y_CurrTrial(ii));
        addpoints(h_eye ,EyeChannel_X_CurrTrial(ii),EyeChannel_Y_CurrTrial(ii));
        
        %{
         if SaveVideo
             
        frame = getframe(gcf);            
        writeVideo(writerObj,frame); 
             
        end
       %}
       
        drawnow limitrate;
  
         clearpoints(h_eye);
           %
        if EraseIndex>EraseLength
           clearpoints(h);
           EraseIndex=0;
        end
        %}
       
        
        %drawnow
        
        EraseIndex=EraseIndex+length(Selector);
        %}
        %{
        if SaveVideo
        frame = getframe(gcf);            
        writeVideo(writerObj,frame); 
        end
        %}
        
        
    end
  
    
clearpoints(h);

 
%{
   if SaveVideo
        frame = getframe(100);            
        writeVideo(writerObj,frame); 
   end 
   %} 
    
    
end
%{
if SaveVideo
close(writerObj);
end
%}




end
