function SceneFractalAssociationCheck(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];



SceneLib=Data.SceneLib;
GoodObjectLib=Data.GoodObjectLib;
BadObjectLib=Data.BadObjectLib;

%{
SceneLib=Data.SceneIndex(:,1);
ObjectLib=Data.TargetIndex(:,:,1);
%}


SceneObjPair=unique([SceneLib(:,3),GoodObjectLib(:,3),BadObjectLib(:,3)],'rows');


SceneData=SceneTuning(Data, StartTrial, EndTrial,0,0);
FractalData=FractalTuning(Data, StartTrial, EndTrial,0,0);

MeanFREachScene=SceneData('AverageFR_EachScene');
SceneName=SceneData('UniqueScene');

GoodFractalFR=FractalData('AverageFR_GoodFractal');
BadFractalFR=FractalData('AverageFR_BadFractal');

FractalName=FractalData('UniqueTarget');


FRScene=MeanFREachScene(:,1);
%FRGoodFractal=GoodFractalFR(:,1);
%FRBadFractal=BadFractalFR(:,1);

FractalNameGood=FractalName(mod(FractalName,floor(FractalName/10))<4);
FractalNameBad=FractalName(mod(FractalName,floor(FractalName/10))>=4);

%First sort scene and fractal firing rates with pair
for i=1:length(SceneName)
    currScene=SceneName(i);
    currGoodObj=SceneObjPair(ismember(SceneObjPair,currScene),2);
    currBadObj=SceneObjPair(ismember(SceneObjPair,currScene),3);
    try
   FRGoodFractal(i,:)=GoodFractalFR(FractalNameGood== currGoodObj,:);
   FRBadFractal(i,:)=BadFractalFR(FractalNameBad== currBadObj,:);
    catch
        keyboard
    end
   
   GoodOrder(i)=find(FractalNameGood== currGoodObj);
   BadOrder(i)=find(FractalNameBad== currBadObj);
   
end

FractalNameGood=FractalNameGood(GoodOrder);
FractalNameBad=FractalNameBad(BadOrder);


%Sort according to the scene
[Val,SortIndex]=sort(FRScene,'descend');

FRScene_Sort=MeanFREachScene(SortIndex,:);
GoodFractalFR_Sort=FRGoodFractal(SortIndex,:);
BadFractalFR_Sort=FRBadFractal(SortIndex,:);

SceneName_Sort=SceneName(SortIndex);
GoodFractalName_Sort=FractalNameGood(SortIndex);
BadFractalName_Sort=FractalNameBad(SortIndex);

if ShowFigureFlag
    FigureStartNum=134;
    FigureIndex=1;
    %Figure 1:Bar plot of scene and fractals
    figtitlestr{FigureIndex}='SceneFractalBarPlot';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
    %First the scenes
    
   subplot(3,1,1);
   g=bar([FRScene_Sort(:,1)],'LineWidth',3);  
  g.FaceColor = 'flat';

ColorChosen= Color(1:length(FRScene_Sort(:,1)),:);
g.CData = ColorChosen;

%ylim([0,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)]); 

hold on

errorbar(FRScene_Sort(:,1),FRScene_Sort(:,2),'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:length(FRScene_Sort(:,1)));
xticklabels(string(SceneName_Sort));
set(gca,'FontSize',15,'FontWeight','Bold')
xlabel('Scene Name')
ylabel('Spikes/s')
box off
%Then the fractals
%Good fractals
 subplot(3,1,2);
   g=bar([GoodFractalFR_Sort(:,1)],'LineWidth',3);  
  g.FaceColor = 'flat';

ColorChosen= Color(1:length(GoodFractalFR_Sort(:,1)),:);
g.CData = ColorChosen;

%ylim([0,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)]); 

hold on

errorbar(GoodFractalFR_Sort(:,1),GoodFractalFR_Sort(:,2),'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:length(GoodFractalFR_Sort(:,1)));
xticklabels(string(GoodFractalName_Sort));
set(gca,'FontSize',15,'FontWeight','Bold')
xlabel('Good Fractal Name')
ylabel('Spikes/s')
box off
  %Bad fractal  
 subplot(3,1,3);
   g=bar([BadFractalFR_Sort(:,1)],'LineWidth',3);  
  g.FaceColor = 'flat';

ColorChosen= Color(1:length(BadFractalFR_Sort(:,1)),:);
g.CData = ColorChosen;

%ylim([0,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)]); 

hold on

errorbar(BadFractalFR_Sort(:,1),BadFractalFR_Sort(:,2),'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:length(BadFractalFR_Sort(:,1)));
xticklabels(string(BadFractalName_Sort));
set(gca,'FontSize',15,'FontWeight','Bold')
xlabel('Bad Fractal Name')
ylabel('Spikes/s')
box off    
    
end
%pause;


end