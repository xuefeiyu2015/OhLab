function [p,rsquare,stats]=anova1_group_stats(x,group,displayflag)
%Function to enable one way anova columnwise
%Xuefei Yu; 02142021
for i=1:size(x,2)
    
   [ p(i),tab,stats{i}]=anova1(x(:,i),group,'off');
   ss=tab{2,2};
   ts=tab{4,2};
   rsquare(i)=ss/ts;
    
end


end
