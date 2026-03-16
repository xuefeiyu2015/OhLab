function  Data_Grouped=SortDataAccordingToGroup(DataBase,SortMarker,NumberOfGroup_Ori,TypeStr);
%Function to sort data according to the markers;
if iscell(NumberOfGroup_Ori)
    GroupSepBin=sort(NumberOfGroup_Ori{1},'descend');
    NumberOfGroup=length(GroupSepBin)+1;
    GroupSepBinUpper=[Inf,GroupSepBin];
    GroupSepBinLower=[GroupSepBin,-Inf];
    
    BoundryFlag=1;
else
    NumberOfGroup=NumberOfGroup_Ori;
    BoundryFlag=0;
end

%Xuefei Yu, 05132020
%Pre: Screen out the NaNs;

SortMarker_NaN=isnan(SortMarker);
DataBase_NaN=DataBase(SortMarker_NaN,:);


SortMarker=SortMarker(~SortMarker_NaN);
DataBase=DataBase(~SortMarker_NaN,:);


%1. Sort Marker
[Marker_Sort, Marker_SortIndex]=sort(SortMarker,TypeStr);
%2. Change the order of the database
Sorted_DataBase=DataBase( Marker_SortIndex,:);

%3.Separate database into groups
TrialNumGroup=ceil(size(DataBase,1)/NumberOfGroup);

LineData=[];
GrouopIndicator=[];

for i=1:NumberOfGroup
    if BoundryFlag
        Data_temp=Sorted_DataBase(Marker_Sort<GroupSepBinUpper(i) & Marker_Sort>=GroupSepBinLower(i),:);
         Marker_temp=Marker_Sort(Marker_Sort<GroupSepBinUpper(i) & Marker_Sort>=GroupSepBinLower(i));
        
    else
   if i~=NumberOfGroup
       Data_temp=Sorted_DataBase((i-1)*TrialNumGroup+1:i*TrialNumGroup,:);
       Marker_temp=Marker_Sort((i-1)*TrialNumGroup+1:i*TrialNumGroup);
       
   else
       Data_temp=Sorted_DataBase((i-1)*TrialNumGroup+1:end,:);
       Marker_temp=Marker_Sort((i-1)*TrialNumGroup+1:end);
   end
    end
   
    Data_Grouped_mean(i,:)=nanmean(Data_temp,1);
     Data_Grouped_sem(i,:)=nanstd(Data_temp,[],1)./sqrt(sum(~isnan(Data_temp)));
    Data_Grouped_Raw{i}=Data_temp;
    DataNum_Grouped(i)=size(Data_temp,1);
    
    MarkerGouped_Raw{i}=Marker_temp;
    MarkerGouped_mean(i)=nanmean(Marker_temp);
    MarkerGouped_sem(i)=nanstd(Marker_temp)/sqrt(sum(~isnan(Marker_temp)));
    %For Stats
    LineData=[LineData;Data_temp];
    GrouopIndicator=[GrouopIndicator;i*ones(size(Data_temp,1),1)];
    
   
end

%Add nans
 Data_Grouped_mean(i+1,:)=nanmean(DataBase_NaN,1);
 Data_Grouped_sem(i+1,:)=nanstd(DataBase_NaN,1)./sqrt(sum(~isnan(DataBase_NaN),1));
 Data_Grouped_Raw{i+1}=DataBase_NaN;
 DataNum_Grouped(i+1)=size(DataBase_NaN,1);
 
Data_Grouped.MarkerMean(i+1)=NaN;
Data_Grouped.MarkerSem(i+1)=NaN;
Data_Grouped.MarkerRaw{i+1}=NaN;


Data_Grouped.mean=Data_Grouped_mean;
Data_Grouped.sem=Data_Grouped_sem;
Data_Grouped.Raw=Data_Grouped_Raw;

Data_Grouped.MarkerMean=MarkerGouped_mean;
Data_Grouped.MarkerSem=MarkerGouped_sem;
Data_Grouped.MarkerRaw=MarkerGouped_Raw;
Data_Grouped.SortIndex=Marker_SortIndex;

Data_Grouped.Num=DataNum_Grouped;
Data_Grouped.LineData=LineData;
Data_Grouped.GrouopIndicator=GrouopIndicator;






end