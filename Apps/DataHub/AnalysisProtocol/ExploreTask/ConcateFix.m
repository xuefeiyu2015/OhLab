function FixTime_Concat=ConcateFix(FixTime);
%Reprocessing fixation time 
%Xuefei Yu 09062020
%Remove the first Fixation if there are more than 2
FixTime_Concat=cellfun(@RemoveFirst,FixTime,'uniform',0);

end
function FT=RemoveFirst(FP_Ori)
Num_F=size(FP_Ori,2);
if Num_F>1
    FT_tmp=FP_Ori(:,2:end);
    FT=[FT_tmp(1,1);FT_tmp(end,end)];
else
    FT=FP_Ori;
end
end