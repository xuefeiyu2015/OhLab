function  Batch_StableValue_Behavior(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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
for i=1:length(FilesName)
  clear OutputData;
  
     load(FilesName{i}{1});
     
     OutputDataTemp=OutputData;
     
     
     Task(i)=OutputData.StableValueBehavior.TaskCode;
     TaskName{i}=OutputData.StableValueBehavior(1).Task;
    
     
     Data =OutputData.StableValueBehavior.DataStamp;

     
     RT_Good{i} = Data('RT_Good')';
     RT_Bad{i} = Data('RT_Bad')';

    RT_Good_Mean(i) = Data('RT_Good_Mean');
    RT_BadMean(i) = Data('RT_Bad_Mean');
     


     
end


RT_Good = cell2mat(RT_Good);
RT_Bad = cell2mat(RT_Bad);

RT_Good_All_Mean = nanmean(RT_Good);
RT_Bad_All_Mean = nanmean(RT_Bad);


RT_Good_All_Sem = nanstd(RT_Good )/sqrt(length(RT_Good ));
RT_Bad_All_Sem = nanstd(RT_Bad)/sqrt(length(RT_Bad));

[h,p]=ttest2(RT_Good,RT_Bad);



if ShowFigureFlag

    figtitlestr{1}='RT_GoodBad';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

b = bar([1,2],[RT_Good_All_Mean,RT_Bad_All_Mean]);
b.FaceColor = 'flat'; % Allow individual face coloring
b.CData(1, :) = [1 0 0]; % Red for the first bar (RGB values)
b.CData(2, :) = [0 0 1]; % Blue for the second bar (RGB values)


hold on
errorbar([1,2],[RT_Good_All_Mean,RT_Bad_All_Mean],[RT_Good_All_Sem,RT_Bad_All_Sem],'sk');
box off
ylabel('Reaction Time(ms)');


end




end


