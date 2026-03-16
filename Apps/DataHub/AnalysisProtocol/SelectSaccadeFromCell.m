function Saccade=SelectSaccadeFromCell(EyeDataXCell,EyeDataYCell,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);
%Updated from SelectSaccade,
%the input x, y are from raw eye trace
%Xuefei Yu 04202020


for i=1:length(EyeDataXCell)
    EyeDataX=EyeDataXCell{i};
    EyeDataY=EyeDataYCell{i};
    if isempty(EyeDataX)
        Saccade(i).NumOfSaccade=0;
         Saccade(i).SaccadeIntervalIndex=NaN;
    
    Saccade(i).Time=NaN;
    Saccade(i).Ecc=NaN;
    
     Saccade(i).X=NaN;
     Saccade(i).Y=NaN;
     Saccade(i).V=NaN;
     Saccade(i).A=NaN;
     Saccade(i).PeakV= NaN;
     Saccade(i).PeakVTime=NaN;
     Saccade(i).SaccadeStartPoint=NaN;
     Saccade(i).SaccadeEndPoint=NaN;
     Saccade(i).SaccadeVector=NaN;
     Saccade(i).SaccadeAngle=NaN;
     Saccade(i).SaccadeAmplitude=NaN;
     Saccade(i).SaccadeStartTime=NaN;
     Saccade(i).SaccadeEndTime=NaN;
     Saccade(i).SaccadeDuration=NaN;
    else
    
    
Ecc=sqrt(EyeDataX.^2+EyeDataY.^2);

%Get velocity from the Eyedata
Ecc_V=diff(Ecc,1,2)/EyeBinWidth*1000;
Ecc_V=[Ecc_V(:,1) ,Ecc_V];

Ecc_Acc=diff(Ecc_V,1,2)/EyeBinWidth*1000;
Ecc_Acc=[Ecc_Acc(:,1),Ecc_Acc];



%First find out the eccentricity more than 3 degree of the base line
%Selection_Ecc=abs(Ecc_Normed)>3;
Selection_V_Low=abs(Ecc_V)>V_Threshold;
%Selection_A_Low=abs(Ecc_Acc)>10000;

%SelectionSaccade_Basic=Selection_Ecc;%.*Selection_V_Low;%.*Selection_A_Low;
SelectionSaccade_Basic=Selection_V_Low;

%figure;

    currEcc=Ecc;
    
    currEcc_X=EyeDataX;
    currEcc_Y=EyeDataY;
    
    currEcc_V=abs(Ecc_V);
    currEcc_A=abs(Ecc_Acc);
    
    currEyeTime=(0:length(EyeDataX)-1)*EyeBinWidth;
    
    
    
    currSelection_Basic=SelectionSaccade_Basic==1;
    
    
    
    SaccadeStartTime=currEyeTime(currSelection_Basic);
    diffStartTime=diff(SaccadeStartTime);
    SeparationPoint=find(diffStartTime>(EyeBinWidth*1.5));
    SeparationPoint=[0,SeparationPoint];
    SaccadeIntervalIndex=[(SeparationPoint(1:end-1)+1)',[SeparationPoint(2:end)]'];
  
    SaccadeIntervalIndex((SaccadeIntervalIndex(:,2)-SaccadeIntervalIndex(:,1))<ContinueBin,:)=[];
    NumberOfSaccade=size(SaccadeIntervalIndex,1);
    
    if NumberOfSaccade>0
    
    %Select out saccade
    SaccadeSelection_Time=repmat(SaccadeStartTime,NumberOfSaccade,1);
 
    SaccadeSelection_Time=SelectMatrixData(SaccadeSelection_Time,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
    
    %Find out the saccade end time
    
    SaccadeStartTime=SaccadeSelection_Time(:,1);
    
    temp=~isnan(SaccadeSelection_Time);
    
    Indices=arrayfun(@(x) find(temp( x,:), 1, 'last'), 1:size(SaccadeSelection_Time, 1));
    SaccadeEndTime=arrayfun(@(x,y) SaccadeSelection_Time(x,y),1:size(SaccadeSelection_Time, 1), Indices)';
    
     %Check saccade interval

     if NumberOfSaccade>1
         InterSaccadeInterval=SaccadeStartTime(2:end)-(SaccadeEndTime(1:end-1));
         
         while (InterSaccadeInterval<ISI_Threshold) & (NumberOfSaccade>1)
             CombinedSaccade=find(InterSaccadeInterval<ISI_Threshold);
             CombineIndex=CombinedSaccade(1);
             
             SaccadeStartTime(CombineIndex+1)=SaccadeStartTime(CombineIndex);
             SaccadeStartTime(CombineIndex)=[];
             SaccadeEndTime(CombineIndex)=[];
             
         %{
             mat_combine=NaN*ones(size(SaccadeSelection_Time));
             mat_combine(CombineIndex,:)=SaccadeSelection_Time(CombineIndex+1,:);
             
             SaccadeSelection_Time=[SaccadeSelection_Time,mat_combine];
             SaccadeSelection_Time(CombineIndex+1,:)=[];
           %}  
             InterSaccadeInterval=SaccadeStartTime(2:end)-(SaccadeEndTime(1:end-1));
             NumberOfSaccade=NumberOfSaccade-1;
             
             %{
             SaccadeIntervalIndex(CombineIndex+1,1)= SaccadeIntervalIndex(CombineIndex,1);
            SaccadeIntervalIndex(CombineIndex,:)=[];
             %}
             
             
         end
         
         
         
         
         
 
     end
     %Select relavent parameters
     
     SaccadeSelection_Time=fillInterval(SaccadeStartTime,SaccadeEndTime,EyeBinWidth);
     SaccadeIntervalIndex=[ceil(SaccadeStartTime/EyeBinWidth)+1,ceil(SaccadeEndTime/EyeBinWidth)+1];
     
     SaccadeSelection_Ecc=repmat(currEcc,NumberOfSaccade,1);
     SaccadeSelection_Ecc=SelectMatrixData(SaccadeSelection_Ecc,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
    
     
     SaccadeSelection_V=repmat(currEcc_V,NumberOfSaccade,1);
     SaccadeSelection_V=SelectMatrixData(SaccadeSelection_V,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     SaccadeSelection_X=repmat(currEcc_X,NumberOfSaccade,1);
     SaccadeSelection_X=SelectMatrixData(SaccadeSelection_X,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     SaccadeSelection_Y=repmat(currEcc_Y,NumberOfSaccade,1);
     SaccadeSelection_Y=SelectMatrixData(SaccadeSelection_Y,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     
     SaccadeSelection_A=repmat(currEcc_A,NumberOfSaccade,1);
     SaccadeSelection_A=SelectMatrixData(SaccadeSelection_A,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     
     %{
     %Select relavent parameters
     
     SaccadeSelection_Ecc=repmat(currEcc(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_Ecc=SelectMatrixData(SaccadeSelection_Ecc,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     
     SaccadeSelection_V=repmat(currEcc_V(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_V=SelectMatrixData(SaccadeSelection_V,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     SaccadeSelection_X=repmat(currEcc_X(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_X=SelectMatrixData(SaccadeSelection_X,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     SaccadeSelection_Y=repmat(currEcc_Y(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_Y=SelectMatrixData(SaccadeSelection_Y,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     
     
     SaccadeSelection_A=repmat(currEcc_A(currSelection_Basic),NumberOfSaccade,1);
     SaccadeSelection_A=SelectMatrixData(SaccadeSelection_A,SaccadeIntervalIndex(:,1),SaccadeIntervalIndex(:,2));
     %}
    %Find out the peak velocity and time
 
     [PeakV,PeakV_Loc] = max(SaccadeSelection_V,[],2);
    
     PeakV_Time=arrayfun(@(x,y) SaccadeSelection_Time(x,y),1:size(SaccadeSelection_Time, 1), PeakV_Loc');
     
     
    %Find out the saccade amplitude, and anlge
   
    SaccadeStartPoint_X=SaccadeSelection_X(:,1);
    SaccadeStartPoint_Y=SaccadeSelection_Y(:,1);
 
    temp=~isnan(SaccadeSelection_X);
    
    Indices=arrayfun(@(x) find(temp( x,:), 1, 'last'), 1:size(SaccadeSelection_X, 1));
    SaccadeEndPoint_X=arrayfun(@(x,y) SaccadeSelection_X(x,y),1:size(SaccadeSelection_X, 1), Indices)';
    SaccadeEndPoint_Y=arrayfun(@(x,y) SaccadeSelection_Y(x,y),1:size(SaccadeSelection_Y, 1), Indices)'; 
   
    SaccadeVector=[SaccadeEndPoint_X,SaccadeEndPoint_Y]-[SaccadeStartPoint_X,SaccadeStartPoint_Y];
  
    %Angle and amplitude
    
   [SaccadeAngle, SaccadeAmplitude]=LocToAngle(SaccadeVector(:,1),SaccadeVector(:,2));
    
   
   
   
   Saccade(i).NumOfSaccade=NumberOfSaccade;
    Saccade(i).SaccadeIntervalIndex=SaccadeIntervalIndex;
    
    Saccade(i).Time=SaccadeSelection_Time;
    Saccade(i).Ecc=SaccadeSelection_Ecc;
    
     Saccade(i).X=SaccadeSelection_X;
     Saccade(i).Y=SaccadeSelection_Y;
     Saccade(i).V=SaccadeSelection_V;
     Saccade(i).A= SaccadeSelection_A;
     Saccade(i).PeakV= PeakV;
     Saccade(i).PeakVTime=PeakV_Time;
     Saccade(i).SaccadeStartPoint=[SaccadeStartPoint_X,SaccadeStartPoint_Y];
     Saccade(i).SaccadeEndPoint=[SaccadeEndPoint_X,SaccadeEndPoint_Y];
     Saccade(i).SaccadeVector=SaccadeVector;
     Saccade(i).SaccadeAngle=SaccadeAngle;
     Saccade(i).SaccadeAmplitude=SaccadeAmplitude;
     Saccade(i).SaccadeStartTime=SaccadeStartTime;
     Saccade(i).SaccadeEndTime=SaccadeEndTime;
     Saccade(i).SaccadeDuration=SaccadeEndTime-SaccadeStartTime;
    else
        Saccade(i).NumOfSaccade=NumberOfSaccade;
         Saccade(i).SaccadeIntervalIndex=NaN;
    
    Saccade(i).Time=NaN;
    Saccade(i).Ecc=NaN;
    
     Saccade(i).X=NaN;
     Saccade(i).Y=NaN;
     Saccade(i).V=NaN;
     Saccade(i).A=NaN;
     Saccade(i).PeakV= NaN;
     Saccade(i).PeakVTime=NaN;
     Saccade(i).SaccadeStartPoint=NaN;
     Saccade(i).SaccadeEndPoint=NaN;
     Saccade(i).SaccadeVector=NaN;
     Saccade(i).SaccadeAngle=NaN;
     Saccade(i).SaccadeAmplitude=NaN;
     Saccade(i).SaccadeStartTime=NaN;
     Saccade(i).SaccadeEndTime=NaN;
     Saccade(i).SaccadeDuration=NaN;
        
    end
     

   
    
    end  
    
end
end
function mat=fillInterval(Start,End,Bin)
MaxCol=max(ceil((End-Start)/Bin)+1);
mat=NaN*ones(size(Start,1),MaxCol);

 for i=1:size(Start,1)
     Number=ceil((End(i,:)-Start(i,:))/Bin)+1;
     mat(i,1:Number)=Start(i,:):Bin:End(i,:);
  
 end





end