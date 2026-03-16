function p=anova1_group(x,group,displayflag)
%Function to enable one way anova columnwise
%Xuefei Yu; 04042020
for i=1:size(x,2)
    
    p(i)=anova1(x(:,i),group,'off');
    
end


end
