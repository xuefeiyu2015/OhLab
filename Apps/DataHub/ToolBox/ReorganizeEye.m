function organized=ReorganizeEye(data);
organized=[];

for i=1:length(data)
    data_curr=data{i};
    try
    maxnumel=max(cellfun(@numel,data_curr));
    
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    catch
    organized{i}=data_curr;
    end
    
end


 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');

end