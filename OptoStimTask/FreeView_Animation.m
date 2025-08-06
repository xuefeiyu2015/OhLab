function FreeView_Animation(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;
TaskType=Data.TaskType;

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
%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);


%Select the eye trace during the 100pre-stim, 100 stim, 100 post stim
%period


%StimStartBin=FindOutTime(EventChannel,EventBin,PRE_STIM);
StimStartBin=FindOutTime(EventChannel,EventBin,ON_STIM )-100;
StimEndBin=FindOutTime(EventChannel,EventBin,PST_STIM );

EyeChannel_X_Stim= SelectEyeInterval(EyeChannel_X_L,StimStartBin,StimEndBin)';
EyeChannel_Y_Stim= SelectEyeInterval(EyeChannel_Y_L,StimStartBin,StimEndBin)';

EyeChannel_X_Stim=ReorganizeEye(EyeChannel_X_Stim);
EyeChannel_Y_Stim=ReorganizeEye(EyeChannel_Y_Stim);

EyeTime_Stim = 0:EyeBinWidth:size(EyeChannel_X_Stim,2)-EyeBinWidth;




MarkerTime = FindOutTime(EventChannel,EventBin,ON_STIM)-StimStartBin;

Interval = [-100,nanmean(StimEndBin-FindOutTime(EventChannel,EventBin,ON_STIM))];
%Interval=[-100,200];
StimDur = [0,nanmean(FindOutTime(EventChannel,EventBin,OFF_STIM)- FindOutTime(EventChannel,EventBin,ON_STIM))];



EyeChannel_X_StimAlign = AlignEyeData(EyeChannel_X_Stim,EyeTime_Stim,MarkerTime,Interval);
EyeChannel_Y_StimAlign = AlignEyeData(EyeChannel_Y_Stim,EyeTime_Stim,MarkerTime,Interval);




%NormX = nanmean(EyeChannel_X_StimAlign(:,1:100),2);
%NormY = nanmean(EyeChannel_Y_StimAlign(:,1:100),2);

%EyeChannel_X_StimAlign_Norm = EyeChannel_X_StimAlign;%-NormX;
%EyeChannel_Y_StimAlign_Norm = EyeChannel_Y_StimAlign;%-NormY;

EyeTime_Stim_Align=Interval(1):EyeBinWidth:Interval(2)+EyeBinWidth;





eyex = reshape(EyeChannel_X_StimAlign',1,numel(EyeChannel_X_StimAlign));

eyey = reshape(EyeChannel_Y_StimAlign',1,numel(EyeChannel_Y_StimAlign));

eyex = eyex(~isnan(eyex));
eyey = eyey(~isnan(eyey));

eyetime_mat =repmat(EyeTime_Stim_Align,size(EyeChannel_X_StimAlign,1),1);
eyetime_line = reshape(eyetime_mat',1,numel(eyetime_mat));



StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
%ControlTrial = ~StimTrial;



trialmarker = repmat(StimTrial,1,size(eyex,2));
trial_line = reshape(trialmarker',1,numel(trialmarker));



if ShowFigureFlag
% Create a figure
f=PrepareFigure(1,[0,0,0],[100,100, 400,400]);
axis([-20 40 -40 40]);
set(gca,'color','k');
hold on;

% Animation loop
for frame = 1:size(eyex,2)% Adjust the number of frames as needed
    %{
    % Update dot positions
    positions = positions + velocities;

    % Wrap dots around the screen
    positions(positions < 0) = 1;
    positions(positions > 1) = 0;
    %}
    
    % Clear the previous frame
    clf
   
    %{
    % Plot the updated dot positions
    scatter(positions(:, 1), positions(:, 2), 10, 'filled');
    
    % Set axis limits
    axis([0 1 0 1]);
    %}
   
    
    % Plot the lines between dots
  %  line(lines(:, [1, 3])', lines(:, [2, 4])');
  dsize = 200;


  timepoint=eyetime_line(frame);
  if timepoint>=0 & timepoint<=100 & trial_line(frame)==1
      eyecolor = 'b';
  else
      eyecolor = 'w';
  end

    scatter(eyex(frame),eyey(frame),dsize,eyecolor,'filled');
    
   axis([-40 40 -40 40]);
   set(gca,'color','k');
   axis off;
    
    % Pause for a short duration to control animation speed
    pause(0.00001);

    
    % Force drawing of the current frame
    drawnow;
   
end %End of the animation





end %End of show figure flag


end %End of the function


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