function [p,h]=ttest_group(x,y,tail)
%Function to enable two sample t test columnwise
%Xuefei Yu; 04042020

p=1;
h=0;

if size(x,2)~=size(y,2)
    
    y=repmat(y,1,size(x,2));
    
end




for i=1:size(x,2)
    [h(i),p(i)]=ttest(x(:,i),y(:,i),'tail',tail);
    
end





end