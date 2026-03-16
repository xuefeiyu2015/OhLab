function Batch_MemorySaccadeLatency(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);

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
ChannelFull=[];
index = 1;
for i=1:length(FilesName)
    for j = 1:length(FilesName{i})
        clear OutputData;
        ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];
  
     load(FilesName{i}{j});
     
     OutputDataTemp=OutputData;
     Task(index)=OutputData.MemorySaccadeLatency.TaskCode;
     TaskName{index}=OutputData.MemorySaccadeLatency.Task;

     FilesNameAll{index}=FilesName_Each{i};
     
     Data =OutputData.MemorySaccadeLatency.DataStamp;

    % PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF');
    % PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF');

    if length(Data('RF_h_v'))>1
        disp('More than 1 RF_hv detected,use the first one instead');
    end
    tmp = Data('RF_h_v');
     RF_h_v(index) =tmp(1);
     tmp = Data('RF_h_m');
     RF_h_m(index)=tmp(1);


     CellType(index) = Data('CellType');
     Latency(index,:) = Data('latencyTime');
     index=index+1;




    end
     
     
end
%PSTH_Time_TgtOn = Data('PSTH_Time_TgtOn');
%PSTH_Time_SacOn = Data('PSTH_Time_SacOn');

VisCells = RF_h_v==1;
L = Latency(VisCells,:);

WhiteW= L(:,1)/1000;
YuW = L(:,2)/1000;


NaNWh = sum(isnan(WhiteW));
NanYu = sum(isnan(YuW));



White=WhiteW(~isnan(WhiteW));
Yu=YuW(~isnan(YuW));

NumTotalW = length(White);
NumTotalY = length(Yu);

NumTotalW_W = length(WhiteW);
NumTotalY_W = length(YuW);


MidWh = median(White);
MidYu = median(Yu);

binWhite= [min(White):0.0001:max(White)+0.0001];

histWhite = histcounts(White,binWhite);

binYu= [min(Yu):0.0001:max(Yu)+0.0001];
histYu = histcounts(Yu,binYu);

cumWhite = cumsum(histWhite)/length(White);
cumYu = cumsum(histYu)/length(Yu);

figure
set(gcf,'color','w');
%h = histogram(White,'Normalization','cumcount','DisplayStyle','stairs');
subplot(2,1,1)
stairs(binWhite(1:end-1),cumWhite,'-k','LineWidth',3);
xlim([0,0.1]);
box off;
grid on;

hold on
plot([0,MidWh],[0.5,0.5],'--r');
plot([MidWh,MidWh],[0,0.5],'-r');

%sprintf('Median=%1.4f',MidWh)
text(MidWh+0.001,0.45,sprintf('Median=%1.4f',MidWh),'color','r');
text(MidWh+0.001,0.35,sprintf('Latency detected for %d of %d neurons',NumTotalW,NumTotalW_W),'color','k');


title('method-white2017');
ylabel('Cumulative proportion of neurons')
xlabel('Time from stimulus onset (s)');

subplot(2,1,2)
stairs(binYu(1:end-1),cumYu,'-k','LineWidth',3);
xlim([0,0.1]);
box off;
grid on;

hold on
plot([0,MidYu],[0.5,0.5],'--r');
plot([MidYu,MidYu],[0,0.5],'-r');

%sprintf('Median=%1.4f',MidWh)
text(MidYu+0.001,0.45,sprintf('Median=%1.4f',MidYu),'color','r');
text(MidYu+0.001,0.35,sprintf('Latency detected for %d of %d neurons',NumTotalY,NumTotalY_W),'color','k');

title('method-yuKatz2023');
ylabel('Cumulative proportion of neurons')
xlabel('Time from stimulus onset (s)');

annotation('textbox',[.05 .7 .3 .3],'String','XF SC during MemorySaccade','LineStyle','none');






end


