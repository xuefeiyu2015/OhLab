function CI=ConfidenceInterval(Data,p)
%To calculate the p% interval for a group of data columnwised
%Data is a matrix, p is in [0,100];
%Xuefei Yu, 03,22,2020
CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x(:)))-1) + mean(x(:),'omitnan'); 

for i=1: size(Data,2)
    x=Data(:,i);
    CI(:,i)=CIFcn(x,p);
    
    %{
    x=Data(~isnan(Data(:,i)),i);
    
    x_Sort=sort(x);
    Num_X=length(x);
    if Num_X<3
        CI(:,i)=[NaN,NaN]';
    else
        Num_High=min(ceil(Num_X*p/100),Num_X);
        Num_Low=max(1,floor(Num_X*(1-p/100)));
        
        CI(:,i)=[x_Sort(Num_Low),x_Sort(Num_High)]';
        
    end
    
    %}
   
   
end




end