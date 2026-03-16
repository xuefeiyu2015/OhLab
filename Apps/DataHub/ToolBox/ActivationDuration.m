function dur = ActivationDuration(bin,resp1,resp2)

dur = 0;
for i = 1:size(resp1,2)
    try
        p(i) = ranksum(resp1(:,i),resp2(:,i));
    catch
        p(i) = NaN;
    end
    
    if p(i) < 0.05
        dur =dur + bin;
    end
    
    
end

end