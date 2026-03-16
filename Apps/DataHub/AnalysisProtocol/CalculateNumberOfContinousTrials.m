function con = CalculateNumberOfContinousTrials(data,errtol);
%Function to calculate the number of contious trials
%error: error tolerance; the maximum error allowed during the continous
%trials(after 5 continous correct trials);
sepindex = find(data==0);
if length(sepindex) > errtol
    diffseq = diff(sepindex);
    %Find out all index that not differ by 1
    conseq = diffseq(diffseq>1)-1;
    startindex_all = sepindex(1:end-1);
    startindex = startindex_all(diffseq>1)+1;
    %Add the last one
    
    lastcon = length(sepindex(end)+1:size(data,2));
    
    if lastcon > 0
        conseq = [conseq,lastcon];
        startindex = [startindex,sepindex(end)+1];
    end
    %Combine the correct segment within errortol
    endindex = conseq + startindex - 1;
    sel = conseq >= 5;
    ori_index = 1:length(conseq);
    if ~isempty(sel)
        
       ori_index_tmp = ori_index(sel);
       %Check if there is any segment to combine
       start_sel = startindex(sel);
       end_sel = endindex(sel);
       com = start_sel(2:end)-end_sel(1:end-1)-1 <= errtol;
       while sum(com)>0
           combine = find(com==1,1,'first');
           ori_index_com =ori_index_tmp(combine);
           
           %Update start_index, end_index
           startindex(ori_index_com+1) =[];
           endindex(ori_index_com) =[];
           conseq(ori_index_com) = conseq(ori_index_com)+conseq(ori_index_com+1);
           conseq(ori_index_com+1)=[];
           
           
           start_sel(combine+1)=[];
           end_sel(combine)=[];
           
           if length(start_sel)<2
               break;
           end
           
           
           
           com = start_sel(2:end)-end_sel(1:end-1)-1 <= errtol;
           
       end
       
        
        
    end
    
    if isempty(conseq)
        conseq = 0;
    end
    con.conseq = conseq;
    con.starttrial = startindex;
    con.endtrial = endindex;
    
else
    con.conseq = sum(data==1);
    con.starttrial = 1;
    con.endtrial = length(1:size(data,2));
    
end



end