function Latency=FindRejectLatency(Dist,Threshold)
%Find out the reject latency
%Xuefei Yu
%05032020

Latency=[];
%{
if length(Threshold)<size(Dist,1)
    Threshold_use=Threshold*ones(size(Dist,1),1);
else
    Threshold_use=Threshold;
end
%}

for i=1:size(Dist,1)
   
  
    dist_curr=Dist(i,:);
    Flag=find(dist_curr>dist_curr(1)+0.5);
    %if sum(Flag)==0
    if isempty(Flag) 
        Latency(i)=NaN;
    else
        Selection=fliplr(dist_curr(1:Flag(1)));
        
        Flag_Update=find( Selection<Threshold);
        if isempty(Flag_Update)
             Flag_Update=find( Selection<dist_curr(1));
            
                 Flag_Update=find( Selection<((dist_curr(1)+ Selection(1))/2));
             
             Latency(i)=Flag(1)-Flag_Update(1)+1;
             
             
            
        
        else
            if Flag(1)-Flag_Update(1)+1<20
                Flag_Update_New=find( Selection<dist_curr(1));
                if isempty(Flag_Update_New)
                    Latency(i)=Flag(1)-Flag_Update(1)+1;
                else
                    Latency(i)=Flag(1)-Flag_Update_New(1)+1;
                end
                
                
            else
                Latency(i)=Flag(1)-Flag_Update(1)+1;
   
            end
        end
        
        
       
 
    end
    
    %Starting Point correction
    if dist_curr(1)>7.5 %Temporally hardwired
       StartPoint=find(dist_curr(1:50)<Threshold);
       if ~isempty(StartPoint)
          Latency(i)=Latency(i)- StartPoint(1);
          disp(sprintf('Starting Point Corrected for %d',StartPoint(1)));
          
       end
    
    end
    %Check abnormalties
    %{
    if Latency(i)>700
        Flag=find(dist_curr>Threshold*1.1);
        Selection=fliplr(dist_curr(1:Flag(1)));
        Flag_Update=find( Selection<Threshold);
        Latency(i)=Flag(1)-Flag_Update(1)+1;
        
        
      %  keyboard
    end
    %}
   
    if Latency(i)<30
        %Check whether the eye came back within 500 bins
        Check_dist= dist_curr(Latency(i):min(Latency(i)+500,length(dist_curr)));
        BackPeriod=find(Check_dist<Threshold);
        if ~isempty(BackPeriod)
           Latency(i)=Latency(i)+BackPeriod(end)-1;
       
        end
        
        
        
    end
    
    
end






end