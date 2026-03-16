function ObjectIndex=ReproduceFromEvent(event,code)

for i=1:size(event,1)
    
   tmp=event(i,ismember(event(i,:),code));
    
 %  number(i)=length(tmp);
   ObjectIndex{i}=tmp;
    
end

end