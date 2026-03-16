function Latency = LatencyWMW(time,resp1,resp2,cri)

Latency = NaN;
for i = 1:length(time)
    try
        p(i) = ranksum(resp1(:,i),resp2(:,i));
    catch
        p(i) = NaN;
    end
    
    %Find out latency
    if i >= cri
        %Get previous data
        check_p = p(i-cri+1 : i);
        if sum(check_p < 0.05) == cri & time(i-cri+1 )>=0
            Latency = time(i-cri+1 );
            break;
        end
        
    end
    
    
end

end