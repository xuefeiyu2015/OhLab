function  RejectReactionAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analysis rejection behavior
%Xuefei Yu 05032020
%Unpack all the data base and task settings
ECode;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
DataPath= Data.Path;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel


FileName=Data.FileName;


%Special request
remove_tmp=[];
SceneRemoval=[];
BlockRemoval=[];

TrialIndexRaw=Data.BlockType(:,1);

%Remove the intermediate scene to focus on the pure scenes if needed
%SceneRemoval=[1054:1057,1060:1063];
if ~isempty(SceneRemoval)
    SceneName_tmp=Data.SceneIndex(Selection,3);
    remove_tmp=ismember(SceneName_tmp,SceneRemoval);
      
end

%Remove the not focus block type to focus on the random blocks if needed
%BlockRemoval=[0,10];
if ~isempty(BlockRemoval)
    Block_tmp=Data.BlockType(Selection,2);
    if ~isempty(remove_tmp)
        remove_tmp=remove_tmp | ismember(Block_tmp,BlockRemoval);
    else
        remove_tmp=ismember(Block_tmp,BlockRemoval);
    end
    
   
end
if ~isempty(remove_tmp)
    Selection=Selection(~remove_tmp);
end
TrialIndex=TrialIndexRaw(Selection);


%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');


%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);
if RecordDate<datestr('04,03,2021') && TaskCode==FORAGING_EXPLORE_TASK_ITL
    fix_tmp=BlockType(:,3)>0;
    BlockType(fix_tmp,2)=BlockType(fix_tmp,2)-10;
    
end


BlockType=BlockType(Selection,2);

AcceptOffset=700;%Time after target onset for acception rounds
TargetWindowThreshold=7.5;





%Correct early block type
%Before 03,18,2020, I only use 0 to indicate fix scene sequence and 1 to
%indicate random scene sequence, both with fixed location sequence
if RecordDate<datestr('03,17,2020')
    
   BlockType(BlockType==0)=10;
    BlockType(BlockType==1)=11;
   
    
end
%BlockType(BlockType(:,3)>0,3)=999;





%BlockType(BlockType==11)=10;
%BlockType(BlockType==1)=0;




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

%Get scene and fractals information
SceneName=Data.SceneIndex(Selection,3);
TargetName=Data.TargetIndex(Selection,2:end,2);
TargetIndex=Data.TargetIndex(Selection,2:end,1);

TargetLoc=Data.TargetLoc(Selection,2:end,:);

TargetNum=Data.TargetNum(Selection,2:end);
TargetSize=Data.TargetSize(Selection,2:end);
TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);
TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);

ChoiceOutcome=Data.ChoiceOutcome(Selection,2:end);

%Select the trials at lease has one case of advance eye in

AdvanceInTargetNum=sum(EventChannel==ADVEYE_IN,2);
AdvanceInTargetTrial=AdvanceInTargetNum>0;

strAdvanceIn=fprintf('There are %d Advance-Eye-In trials\n',sum(AdvanceInTargetTrial));


%Select Rejection time period: from target onset to the next fixation onset
 
% TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,ADVEYE_IN);
 TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
  FixOnTime=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD );
%FixOnTime=FindOutTime(EventChannel,EventTimeChannel,REJECT);
 
 %screen out the searching trials
 TargetOnTime=TargetOnTime(AdvanceInTargetTrial,:);
 FixOnTime=FixOnTime(AdvanceInTargetTrial,:);
 
 
 FixOnTime=FixOnTime(:,2:end);%Remove the first fixation onset
 FixOnTime=[FixOnTime,NaN*ones(size(FixOnTime,1),1)];
 
 
 SceneName=SceneName(AdvanceInTargetTrial,:);
 
 EyeChannel_X_L=EyeChannel_X_L(AdvanceInTargetTrial);
 EyeChannel_Y_L=EyeChannel_Y_L(AdvanceInTargetTrial);
 
 TargetLoc=TargetLoc(AdvanceInTargetTrial,:,:);
 TargetPolarAngle=TargetPolarAngle(AdvanceInTargetTrial,:);
 
 TargetIndex=TargetIndex(AdvanceInTargetTrial,:,:);
 TargetNum=TargetNum(AdvanceInTargetTrial,:);
 
 BlockType=BlockType(AdvanceInTargetTrial);
 
 
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};


%Select the eye movements between the target onset and the next fixation onset
for i=1:size(TargetOnTime,2)
    
    StartTimeBin=floor((TargetOnTime(:,i))/EyeBinWidth);   
    
    EndTimeBin=ceil(FixOnTime(:,i)/EyeBinWidth);
    
    temp=floor((~isnan(StartTimeBin))&isnan(EndTimeBin));
    temp(isnan(StartTimeBin))=NaN*temp(isnan(StartTimeBin));
    
    AcceptionFlag(:,i)=temp;
    
    EndTimeBin(isnan(EndTimeBin))=StartTimeBin(isnan(EndTimeBin))+AcceptOffset/EyeBinWidth;
    
    X_temp=SelectEyeInterval(EyeChannel_X_L,StartTimeBin,EndTimeBin);
    Y_temp=SelectEyeInterval(EyeChannel_Y_L,StartTimeBin,EndTimeBin);
    
    EyeChannel_X_RA_Period=X_temp{1};
    EyeChannel_Y_RA_Period=Y_temp{1};
    
    %To mat
    maxnumel=max(cellfun(@numel,EyeChannel_X_RA_Period));
    
    EyeChannel_X_RA_Period=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-numel(x))],EyeChannel_X_RA_Period,'uniform',0));
    EyeChannel_Y_RA_Period=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-numel(x))],EyeChannel_Y_RA_Period,'uniform',0));
    
   
   %Measure the distance
   TargetLoc_curr=squeeze(TargetLoc(:,i,:));
   TargetType=squeeze(TargetIndex(:,i));
   
   GoodObjectFlag=mod(TargetType,2);
   BadObjectFlag=floor(mod(TargetType,2)==0);
   BadObjectFlag(isnan(TargetType))=NaN;
   
   
   DistanceToTarget=sqrt((EyeChannel_X_RA_Period-TargetLoc_curr(:,1)).^2+(EyeChannel_Y_RA_Period-TargetLoc_curr(:,2)).^2);
   
   DistanceToTarget_Accept{i}=DistanceToTarget(temp==1,:);
   DistanceToTarget_Reject{i}=DistanceToTarget(temp==0,:);
   
  
   DistanceToTarget_GoodAccept{i}=DistanceToTarget((temp==1)&(GoodObjectFlag==1),:);
   DistanceToTarget_GoodReject{i}=DistanceToTarget((temp==0)&(GoodObjectFlag==1),:);
   
   DistanceToTarget_BadAccept{i}=DistanceToTarget((temp==1)&(BadObjectFlag==1),:);
   DistanceToTarget_BadReject{i}=DistanceToTarget((temp==0)&(BadObjectFlag==1),:);
   
   
   %Measure reject latency
  % RejectLatency=((FindRejectLatency(DistanceToTarget,DistanceToTarget(:,1))-1)*EyeBinWidth)';%Transform bin to time

    RejectLatency=(((FindRejectLatency(DistanceToTarget,TargetWindowThreshold)-1)*EyeBinWidth))';%Transform bin to time
    
    
    RejectLatency_Accept=RejectLatency(temp==1);
    RejectLatency_Reject=RejectLatency(temp==0);
    
      RejectLatency_BadReject{i}=RejectLatency(temp==0&BadObjectFlag==1);
     RejectLatency_GoodReject{i}=RejectLatency(temp==0&GoodObjectFlag==1);
    
    
 %Store this for example plots  
 EyeChannel_X_RA_Period_BadReject{i}=EyeChannel_X_RA_Period(temp==0&BadObjectFlag==1,:);
 EyeChannel_Y_RA_Period_BadReject{i}=EyeChannel_Y_RA_Period(temp==0&BadObjectFlag==1,:);
   
   
 EyeChannel_X_RA_Period_GoodAccept{i}=EyeChannel_X_RA_Period(temp==1&GoodObjectFlag==1,:);
 EyeChannel_Y_RA_Period_GoodAccept{i}=EyeChannel_Y_RA_Period(temp==1&GoodObjectFlag==1,:);
 
  RejectLatency_All(:,i)= RejectLatency;
    
   
   
    
end

%For histograms
RejectLatency_BadReject_All=cell2mat(cellfun(@(x) x',RejectLatency_BadReject,'uniform',0));
RejectLatency_GoodReject_All=cell2mat(cellfun(@(x) x',RejectLatency_GoodReject,'uniform',0));
Bin=10;
Edges=[0:Bin:600];
BadReject_LatencyCounts = histcounts(RejectLatency_BadReject_All,Edges);
GoodReject_LatencyCounts = histcounts(RejectLatency_GoodReject_All,Edges);

Mean_Reject_LatencyAll=mean([RejectLatency_BadReject_All,RejectLatency_GoodReject_All]);
TotalNumberOfReject=numel(RejectLatency_BadReject_All)+numel(RejectLatency_GoodReject_All);

SD_Reject_LatencyAll=std([RejectLatency_BadReject_All,RejectLatency_GoodReject_All]);%./sqrt(TotalNumberOfReject);


fprintf('Total Number of Rejections: %d',TotalNumberOfReject);



%Choose one  example eye trace and distance for accept and reject

%Set criteria to find out a sample 
AllTargetLoc_ForExample=squeeze(TargetLoc(:,1,:));

GoodObjectFlag=mod(TargetIndex,2);
BadObjectFlag=mod(TargetIndex,2)==0;


%AllTargetLoc_GoodAccept=AllTargetLoc_ForExample((AcceptionFlag(:,1)==1) & (TargetNum==1),:);

AllTargetLoc_GoodAccept=AllTargetLoc_ForExample((AcceptionFlag(:,1)==1)&(GoodObjectFlag(:,1)==1),:);

%AllTargetLoc_BadReject=AllTargetLoc_ForExample((AcceptionFlag(:,1)==0) & (TargetNum>1),:);
AllTargetLoc_BadReject=AllTargetLoc_ForExample((AcceptionFlag(:,1)==0)&(GoodObjectFlag(:,1)==0),:);

AcceptTrialSample=15;
ExampleTargetLoc=AllTargetLoc_GoodAccept(AcceptTrialSample,:);
RejectTrialSampleLib=find(sum(AllTargetLoc_BadReject==ExampleTargetLoc,2)==2);

RejectTrialSample=RejectTrialSampleLib(min(10,size(RejectTrialSampleLib,1)));

Example_Accept_X=EyeChannel_X_RA_Period_GoodAccept{1}(AcceptTrialSample,:);
Example_Accept_Y=EyeChannel_Y_RA_Period_GoodAccept{1}(AcceptTrialSample,:);

Example_Reject_X=EyeChannel_X_RA_Period_BadReject{1}(RejectTrialSample,:);
Example_Reject_Y=EyeChannel_Y_RA_Period_BadReject{1}(RejectTrialSample,:);

ExampleDistanceToTarget_Reject=DistanceToTarget_BadReject{1}(RejectTrialSample,:);
ExampleDistanceToTarget_Accept=DistanceToTarget_GoodAccept{1}(AcceptTrialSample,:);

%Reject Latency For the example
%RejectLatency_Example_temp=RejectLatency_All((AcceptionFlag(:,1)==0) & (TargetNum>1),:);
RejectLatency_Example_temp=RejectLatency_All((AcceptionFlag(:,1)==0) &(GoodObjectFlag(:,1)==0),:);

RejectLatency_Example=RejectLatency_Example_temp(RejectTrialSample,1);
fprintf('Example Reject Latency: %3.2f (ms)',RejectLatency_Example);

%Analysis reject latency for bad objects as a function of the fractal
%order
temp_trans=AcceptionFlag;
temp_trans(AcceptionFlag==1)=NaN;
temp_trans(AcceptionFlag==0)=1;
RejectLatency_Bad_EachOrder= RejectLatency_All.*temp_trans;
%RejectLatency_BadReject_Mat=

Num_BadReject_EachOrder=sum(~isnan(RejectLatency_Bad_EachOrder),1);
RejectLatency_Bad_EachOrder_mean=nanmean(RejectLatency_Bad_EachOrder,1);

%First measure the rejection speed among blocks 
Diff_First2Last=[];
RejectionNum_ForCompare=[];

FirstRejection_All=[];
LastRejection_All=[];

for i=1:NumberOfBlocks
    CurrentBlock=UniqueBlockType(i);
    BlockOrder(i)=find(AssociateBlock==CurrentBlock);
    BadRejectLatency_EachOrder_CurrBlock_mean(i,:)=nanmean(RejectLatency_Bad_EachOrder(BlockType==CurrentBlock,:),1);
    BadRejectLatency_EachOrder_CurrBlock_sem(i,:)=nanstd(RejectLatency_Bad_EachOrder(BlockType==CurrentBlock,:),[],1)./sqrt(sum(isnan(RejectLatency_Bad_EachOrder(BlockType==CurrentBlock,:)),1));
   
    
    %Histogram for each blocks
    
    RejectLatency_BadReject=RejectLatency_Bad_EachOrder(BlockType==CurrentBlock,:);
    RejectLatency_BadReject=reshape(RejectLatency_BadReject,1,numel(RejectLatency_BadReject));
    RejectLatency_BadReject=RejectLatency_BadReject(~isnan(RejectLatency_BadReject));
    
    BadReject_LatencyCounts_Block(i,:) = histcounts(RejectLatency_BadReject,Edges);
    BadReject_LatencyCounts_mean_Block(i)=nanmean(RejectLatency_BadReject);
    BadReject_LatencyCounts_sem_Block(i)=nanstd(RejectLatency_BadReject)/sqrt(sum(~isnan(RejectLatency_BadReject)));
   
    BadReject_LatencyCompare{i}=RejectLatency_BadReject;



   
    
    %To push it further
    Diff_BadReject_Order=diff(RejectLatency_Bad_EachOrder(BlockType==CurrentBlock,:),1,2);
    
    Diff_BadRejectLatency_EachOrder_CurrBlock_mean(i,:)=nanmean(Diff_BadReject_Order,1);
     Diff_BadRejectLatency_EachOrder_CurrBlock_sem(i,:)=nanstd(Diff_BadReject_Order,1)./sqrt(sum(isnan(Diff_BadReject_Order),1));
     
     %Compare the first rejection and the last rejection
     RejectionLib=RejectLatency_Bad_EachOrder(BlockType==CurrentBlock,:);
     RejectionNumLib=TargetNum(BlockType==CurrentBlock,:)-1;
     %Only use trials with >1 cases of rejections
     FirstRejection=RejectionLib(RejectionNumLib>1,1);
     LastRejection=RejectionLib(RejectionNumLib>1,[RejectionNumLib(RejectionNumLib>1):RejectionNumLib(RejectionNumLib>1)]);
     
     FirstRejection_mean(i)=nanmean(FirstRejection);
     LastRejection_mean(i)=nanmean(LastRejection);
     
     FirstRejectionBlock{i}=FirstRejection;
     LastRejectionBlock{i}=LastRejection;
     
     
     
     RejectionNum_ForCompare=[RejectionNum_ForCompare;RejectionNumLib(RejectionNumLib>1)];
     
     
     Diff_First2Last=[Diff_First2Last;[LastRejection- FirstRejection]];
     
     [h,p_First2LastReject(i)]=ttest(FirstRejection,LastRejection);
     
     FirstRejection_All=[FirstRejection_All;FirstRejection];
     LastRejection_All=[LastRejection_All;LastRejection];
     
     
end

%Overall significance

[h,p_FirstRLastR_All]=ttest(FirstRejection_All,LastRejection_All);

PropFirstLarger=sum(LastRejection_All<FirstRejection_All)/length(LastRejection_All);
PropLastLarger=sum(LastRejection_All>=FirstRejection_All)/length(LastRejection_All);

%Each fractal order
%RejectLatency_Bad_EachOrder_mean=RejectLatency_Bad_EachOrder

%figure
%plot(RejectionNum_ForCompare,Diff_First2Last,'or')


%plot(RejectLatency_Bad_EachOrder(ismember(BlockType(:,1),[10,11]),:)','-')

%Overall significance between blocks
if length(BadReject_LatencyCompare)>1
    p_SigLatencyDiff=ranksum(BadReject_LatencyCompare{1},BadReject_LatencyCompare{2});
    
    
end


if ShowFigureFlag
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %Figure 1 :Plot the distribution of Reject Latency for both good object
    %and the bad object
titlestr{1}='RejectLatencyDistribution';
fig{1}=PrepareFigure(131,'w',[50,100, 800,400],'Name',titlestr{1});

b=bar([(Edges(1:end-1)+Bin/2)',(Edges(1:end-1)+Bin/2)'],[GoodReject_LatencyCounts',BadReject_LatencyCounts'],'stacked');
b(1).FaceColor='#DC2648';
b(2).FaceColor='#648FFF';
b(1).EdgeColor='#DC2648';
b(2).EdgeColor='#648FFF';
%Add mean
hold on
yrange=max(BadReject_LatencyCounts+GoodReject_LatencyCounts);
p=plot(Mean_Reject_LatencyAll,yrange*1.1,'V','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');

legend('RejectGood','RejectBad')
textstring=sprintf('Mean: %2.1f',Mean_Reject_LatencyAll);
text(Mean_Reject_LatencyAll*0.8,yrange*1.2,textstring,'FontSize',25,'FontWeight','Bold');
box off
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
xlabel('Rejection Latency (ms)','FontSize',25,'FontWeight','Bold');
ylabel('Number of cases','FontSize',25,'FontWeight','Bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 2: Plot the example eye trances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% and distance dynamics for a few cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


titlestr{2}='Example eye traces and distance for rejection';
fig{2}=PrepareFigure(132,'w',[70,100, 800,500],'Name',titlestr{2});
subplot(1,2,1)
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
                
              
                t=cplot(r,ExampleTargetLoc(1),ExampleTargetLoc(2));
                t.FaceColor=ColorFractal;
                t.FaceAlpha=0.1;
                t.EdgeColor=ColorFractal;
                t.EdgeAlpha=0.1;
                
               %Plot fixation location
         currFixLocX=0;
         currFixLocY=0;
         
         p=plot(currFixLocX,currFixLocY,'sg','MarkerSize',15, 'MarkerFaceColor', 'g');
         %Plot the target center
         plot(ExampleTargetLoc(1),ExampleTargetLoc(2),'xg','MarkerSize',15, 'MarkerEdgeColor', ColorFractal,'LineWidth',3);
         
           
%Acception
plot(Example_Accept_X',Example_Accept_Y','-r','LineWidth',3);
%Emphasize starting point and the end point
plot(Example_Accept_X(1)',Example_Accept_Y(1)','.r','MarkerSize',20,'LineWidth',3);
%plot(Example_Accept_X(end)',Example_Accept_Y(end)','r','Marker','+','MarkerSize',25,'LineWidth',3);
%Rejection
plot(Example_Reject_X',Example_Reject_Y','-b','LineWidth',3);
plot(Example_Reject_X(1)',Example_Reject_Y(1)','.b','MarkerSize',20,'LineWidth',3);



axis equal
xlim([-25,25])
ylim([-25,25])
axis off

%Plot distance
subplot(1,2,2)
%figure
plot(ExampleDistanceToTarget_Reject,'-b','LineWidth',3);
hold on
plot(ExampleDistanceToTarget_Reject(1),'.b','MarkerSize',20,'LineWidth',3);
plot(ExampleDistanceToTarget_Accept,'-r','LineWidth',3);
plot(ExampleDistanceToTarget_Accept(1),'.r','MarkerSize',20,'LineWidth',3);
%ylim([0,20]);
yrange=ylim;
xrange=xlim;
plot([RejectLatency_Example,RejectLatency_Example],[yrange(1),yrange(2)],'--k','LineWidth',3);
plot([xrange(1),xrange(2)],[7.5,7.5],'--k','LineWidth',3);


set(gca,'LineWidth',3,'FontWeight','Bold','FontSize',20);
box off
xlabel('Time from Target Onset(ms)','FontSize',25,'FontWeight','Bold')
ylabel('Distance from the target center','FontSize',25,'FontWeight','Bold')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compare the first rejection and the last rejection
titlestr{3}='FristRejectionVSLastRejection';
fig{3}=PrepareFigure(133,'w',[90,100, 800,800],'Name',titlestr{3});
axes('Position',[0.15,0.1,0.6,0.8]);
plot([FirstRejection_All';LastRejection_All'],'Marker','o','LineStyle','-','Color','#B5B5B7');
hold on
b=bar([1,2],[nanmean(FirstRejection_All),nanmean(LastRejection_All)]);
b.FaceColor = 'none';
b.EdgeColor='k';
b.LineWidth=3;
b.BarWidth=0.5;
hold on
errorbar([1,2],[nanmean(FirstRejection_All),nanmean(LastRejection_All)],[nanstd(FirstRejection_All)/sqrt(length(FirstRejection_All)),nanstd(LastRejection_All)/sqrt(length(LastRejection_All))],'.k','LineWidth',3,'MarkerSize',10)
xlim([0.5,2.5]);
xticks([1,2]);
xticklabels({'FirstRej','LastRej'});

%Add significance
if p_FirstRLastR_All<0.05
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on

yrange=ylim;
plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_FirstRLastR_All<0.05 & p_FirstRLastR_All>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.2, '*k','MarkerSize',10)
elseif p_FirstRLastR_All<0.01 & p_FirstRLastR_Alle>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.2, '*k','MarkerSize',10)
        
elseif p_FirstRLastR_All<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.2, '*k','MarkerSize',10)
    
end
end

set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
ylabel('Reject Latency (ms)','FontSize',30);
title('Comparison in total','FontSize',30);
box off

%Internal pannel
axes('Position',[0.7,0.5,0.23,0.25]);
b2=bar([1,2],[PropFirstLarger,PropLastLarger]*100);
b2.FaceColor = '#999CA4';
b2.EdgeColor='#999CA4';
b2.LineWidth=3;
b2.BarWidth=0.5;

xlim([0.5,2.5]);
xticks([1,2]);
xticklabels({'Decre','Incre'});
ylim([0,100]);
ylabel('Prop of trials');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold')
box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compare the first rejection and the
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%last rejection for each blocks
titlestr{4}='FristRejectionVSLastRejection_Block';
fig{4}=PrepareFigure(134,'w',[110,100, 800,600],'Name',titlestr{4});
for i=1:NumberOfBlocks
     CurrentBlock=UniqueBlockType(i);
   
    
    subplot(ceil(NumberOfBlocks/2),2,i);
    plot([FirstRejectionBlock{i}';LastRejectionBlock{i}'],'Marker','o','LineStyle','-','Color','#B5B5B7');
hold on
b=bar([1,2],[nanmean(FirstRejectionBlock{i}),nanmean(LastRejectionBlock{i})]);
b.FaceColor = 'none';
b.EdgeColor='k';
b.LineWidth=3;
b.BarWidth=0.5;
hold on
errorbar([1,2],[nanmean(FirstRejectionBlock{i}),nanmean(LastRejectionBlock{i})],[nanstd(FirstRejectionBlock{i})/sqrt(length(FirstRejectionBlock{i})),nanstd(LastRejectionBlock{i})/sqrt(length(LastRejectionBlock{i}))],'.k','LineWidth',3,'MarkerSize',10)
xlim([0.5,2.5]);
xticks([1,2]);
xticklabels({'FirstRej','LastRej'});

%Add significance
if p_First2LastReject(i)<0.05
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on

yrange=ylim;
plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_First2LastReject(i)<0.05 & p_First2LastReject(i)>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.2, '*k','MarkerSize',10)
elseif p_First2LastReject(i)<0.01 & p_First2LastReject(i)>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.2, '*k','MarkerSize',10)
        
elseif p_First2LastReject(i)<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.2, '*k','MarkerSize',10)
    
end
end

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
ylabel('Reject Latency(ms)');
title(Blockstr{BlockOrder(i)});
box off
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compare the Rejection speed for each
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr{5}='RejectionSpeed_EachFractal';
fig{5}=PrepareFigure(135,'w',[130,100, 800,600],'Name',titlestr{5});
%RejectLatency_Bad_EachOrder

plot(RejectLatency_Bad_EachOrder','Marker','o','LineStyle','-','Color','#B5B5B7','LineWidth',2);
hold on
b=bar(1:size(RejectLatency_Bad_EachOrder,2),nanmean(RejectLatency_Bad_EachOrder,1));
b.FaceColor = 'none';
b.EdgeColor='k';
b.LineWidth=3;
b.BarWidth=0.5;
hold on
errorbar([1:size(RejectLatency_Bad_EachOrder,2)],nanmean(RejectLatency_Bad_EachOrder,1),nanstd(RejectLatency_Bad_EachOrder)./sqrt(sum(isnan(RejectLatency_Bad_EachOrder))),'.k','LineWidth',3,'MarkerSize',10)
xlim([0.5,size(RejectLatency_Bad_EachOrder,2)+0.5]);
%xticks([1,2]);
%xticklabels({'FirstRej','LastRej'});

%{
%Add significance
if p_FirstRLastR_All<0.05
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on

yrange=ylim;
plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_FirstRLastR_All<0.05 & p_FirstRLastR_All>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.2, '*k','MarkerSize',10)
elseif p_FirstRLastR_All<0.01 & p_FirstRLastR_Alle>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.2, '*k','MarkerSize',10)
        
elseif p_FirstRLastR_All<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.2, '*k','MarkerSize',10)
    
end
end
%}

set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
xlabel('Rejection Order','FontSize',30,'FontWeight','Bold');
ylabel('Reject Latency (ms)','FontSize',30,'FontWeight','Bold');
title('Comparison of Reject Latency for each order','FontSize',30,'FontWeight','Bold');
box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    % :Plot the distribution of Reject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Latency for each block
    
titlestr{6}='RejectLatencyDistribution';
fig{6}=PrepareFigure(136,'w',[150,100, 1200,600],'Name',titlestr{6});


for i=1:NumberOfBlocks
    subplot(ceil(NumberOfBlocks/2),2,i);

b=bar((Edges(1:end-1)+Bin/2),BadReject_LatencyCounts_Block(i,:));
b.FaceColor=ColorForEachMode(BlockOrder(i),:);
b.EdgeColor=b.FaceColor;

%Add mean
hold on
yrange=max(BadReject_LatencyCounts_Block(i,:),[],'all');
p=plot(BadReject_LatencyCounts_mean_Block(i),yrange*1.1,'V','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');


textstring=sprintf('Mean: %2.1f',BadReject_LatencyCounts_mean_Block(i));
text(BadReject_LatencyCounts_mean_Block(i)*0.8,yrange*1.2,textstring,'FontSize',30,'FontWeight','Bold');
ylim([0,yrange*1.4]);
box off
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');

if mod(NumberOfBlocks,2)==0
    if i>NumberOfBlocks-2
        xlabel('Rejection Latency (ms)','FontSize',30,'FontWeight','Bold');
       
    end
else
    if i>NumberOfBlocks-1
         xlabel('Rejection Latency (ms)','FontSize',30,'FontWeight','Bold');
       
    end
    
end
if mod(i,2)==1
     ylabel('Number of cases','FontSize',30,'FontWeight','Bold');
    
end




title(Blockstr{BlockOrder(i)},'FontSize',30,'FontWeight','Bold');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %Compare the overall reject latency

titlestr{7}='RejectLatencyComparisonbetweenBlocks';
fig{7}=PrepareFigure(137,'w',[170,100, 600,600],'Name',titlestr{7}); 

b=bar([1:length(BadReject_LatencyCounts_mean_Block)],BadReject_LatencyCounts_mean_Block);

b.FaceColor='flat';
b.EdgeColor=b.FaceColor;
b.CData=ColorForEachMode(BlockOrder,:);
hold on
errorbar([1:length(BadReject_LatencyCounts_mean_Block)],BadReject_LatencyCounts_mean_Block,BadReject_LatencyCounts_sem_Block,'sk','LineWidth',3);

%Add significance

%p_SigLatencyDiff
  box off
 set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
    ylabel('Mean & SEM','FontSize',25,'FontWeight','Bold');
   
xticklabels([Blockstr{BlockOrder}]);
yrange=ylim;
%Add significance level
p=p_SigLatencyDiff;
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

end 
     



end