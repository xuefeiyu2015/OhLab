function DataGrouped=PoloarHistGroup(Angle,Data,Edge)
%Function to bin on the polar coordinates
%Xuefei Yu
%Angles: in degree
DataGrouped=[];
DataLine=[];
GroupLine=[];

%Columnwise;
for j=1:length(Edge)
    if j==length(Edge)
          Interval=[Edge(j),Edge(1)];
          SelectAngle=(Angle>=Interval(1) & Angle<359.999)|(Angle>=0 & Angle<Interval(2));
          Center(j)=(Edge(1)+(Edge(j)-360))/2;
            
    else
           Interval=[Edge(j),Edge(j+1)];
           SelectAngle=Angle>=Interval(1) & Angle<Interval(2);
          Center(j)=mean(Interval);
    end
    if sum(SelectAngle)<3
        DataInBin=NaN;
        DataInBin_mean(j,:)=NaN*ones(1,size(Data,2));
        DataInBin_sem(j,:)=NaN*ones(1,size(Data,2));
        DataInBin_Num(j)=sum(SelectAngle);
        DataInBin_Raw{j}=NaN;
       % Center(j)=NaN;
        
    else
    
         DataInBin=Data(SelectAngle,:);
         DataInBin_mean(j,:)=nanmean(DataInBin,1);
         DataInBin_sem(j,:)=nanstd(DataInBin,[],1)./sqrt(sum(~isnan(DataInBin),1));
       
         
         DataInBin_Num(j)=sum(SelectAngle);
          DataInBin_Raw{j}=DataInBin;
         
         %For Stats
         DataLine=[DataLine;DataInBin];
         GroupLine=[GroupLine;repmat(Center(j),DataInBin_Num(j),1)];
        
        
     end  
end
    
    
%Stats via ANOVA

  try  
 ANOVA_p=anova1_group(DataLine,GroupLine,'off');
  catch
      ANOVA_p=NaN;
  end
 
 %Prefer direction via Vector Summation
 for ii=1:size(DataInBin_mean,2)
     try
 [PrefVector(ii), ~, PrefVector_Amp(ii)] = vectorsumAngle(DataInBin_mean(:,ii),  Center, 0*ones(length(Center),1));
     catch
         PrefVector(ii)=NaN;
         PrefVector_Amp(ii)=NaN;
     end
 
 end
 
 %Prefer direction based on condition numbers
 
 
     try
 [PrefVectorNum, ~, PrefVector_AmpNum] = vectorsumAngle(DataInBin_Num,  Center, 0*ones(length(Center),1));
     catch
         PrefVectorNum(ii)=NaN;
         PrefVector_AmpNum(ii)=NaN;
     end
 
 
  
 

DataGrouped.DataMean=DataInBin_mean;
DataGrouped.DataSEM=DataInBin_sem;
DataGrouped.DataNum=DataInBin_Num;
DataGrouped.Center=Center;

DataGrouped.DataLine=DataLine;
DataGrouped.GroupLine=GroupLine;
DataGrouped.DataInBin_Raw=DataInBin_Raw;
DataGrouped.ANOVA_p=ANOVA_p;
DataGrouped.PrefVector=PrefVector;
DataGrouped.PrefVector_Amp=PrefVector_Amp;

DataGrouped.PrefVectorNum=PrefVectorNum;
DataGrouped.PrefVector_AmpNum=PrefVector_AmpNum;









%Group=PoloarHistGroup(SaccadeAngle_First,[SpikeRate_PreSac',SpikeRate_DurSac',SpikeRate_PostSac],SeparationPoint)




end