function AlignedSpike=SpikeAlignTool(AlignedTime,SpikeTime_Ori,StepSize,BinWidth);
%Function to align spikes according to different time marker 
%Xuefei Yu 06072020


SpikeTime_Start=SpikeTime_Ori(:,end-1);
SpikeTime_End=SpikeTime_Ori(:,end);

if sum(~(isnan(SpikeTime_Start)))~=size(SpikeTime_Start,1) ||sum(~(isnan(SpikeTime_End)))~=size(SpikeTime_End,1) 
    SpikeTime_Start=0*ones(size(SpikeTime_Start,1),1);
    SpikeTime_End=9999*ones(size(SpikeTime_End,1),1);
    disp('No start/end information, use inf instead');
end


SpikeTime_Aligned_Whole=[];
SpikeCount_Aligned_Whole=[];

for i=1:size(AlignedTime,2)
    SpikeTime_Aligned=SpikeTime_Ori-AlignedTime(:,i);
    SelectIntervalTime=[SpikeTime_Start-AlignedTime(:,i), SpikeTime_End-AlignedTime(:,i)];  
    if sum(~isnan(AlignedTime(:,i)))==0
        SpikeFiringRateAligned{i}=NaN;
        TimeSequence_Aligned{i}=NaN;
        SpikeCountAligned{i}=NaN;
    else
    %Using a gaussian smooth window to calculate firing rate
   
        [SpikeFiringRateAligned{i},TimeSequence_Aligned{i},SpikeCountAligned{i}]=CalculateSpikeFiringRate(SpikeTime_Aligned,StepSize,BinWidth,SelectIntervalTime);
   
    end
        %PSTH_Aligned_Whole=[PSTH_Aligned_Whole;SpikeFiringRateAligned{i}];
      %  PSTH_Aligned_Time_Whole=[PSTH_Aligned_Time_Whole;TimeSequence_Aligned{i}];
      SpikeTime_Aligned_Whole=[SpikeTime_Aligned_Whole;SpikeTime_Aligned];  
     
end


max_length=max(cellfun(@(x) size(x,2), SpikeFiringRateAligned));

PSTH_Aligned=cell2mat(cellfun(@(x)[x,nan*ones(size(x,1),max_length-size(x,2))],SpikeFiringRateAligned,'uniform',0)');
PSTH_Time_Aligned=cell2mat(cellfun(@(x)[x,nan*ones(size(x,1),max_length-size(x,2))],TimeSequence_Aligned,'uniform',0)');

PSTH_SpikeCountAligned=cell2mat(cellfun(@(x)[x,nan*ones(size(x,1),max_length-size(x,2))],SpikeCountAligned,'uniform',0)');

AlignedSpike.SpikeTime_Aligned_Whole=SpikeTime_Aligned_Whole;
AlignedSpike.PSTH_Aligned=PSTH_Aligned;
AlignedSpike.PSTH_Time_Aligned=PSTH_Time_Aligned;
AlignedSpike.PSTH_SpikeCount_Aligned=PSTH_SpikeCountAligned;


end