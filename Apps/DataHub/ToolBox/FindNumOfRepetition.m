function [temp_loc,NumInRep,RepData]=FindNumOfRepetition(data);
%Return the numer of conditions in each repetition
%Xuefei Yu; 04162020
unique_data=unique(data,'stable');

MaxRep=ceil(length(data)/length(unique_data));

data_full=[data;NaN*(ones(MaxRep*length(unique_data)-length(data),1))];

data_reshape=reshape(data_full,length(unique_data),MaxRep);
RepData=NaN*ones(size(data_reshape));

for i=1:MaxRep
    curr_data=data_reshape(:,i);
    CheckRep=countmember(curr_data,unique_data);
    if sum(CheckRep>1)
        loc=find(CheckRep>1);
        RepData(:,i)=curr_data;
        RepData(:,loc:end)=NaN;
        
        temp=[curr_data(:,loc:end)',reshape(1,numel(data_reshape(:,[i+1:end])))];
        MaxRep_temp=ceil(length(temp)/length(unique_data));
        
        data_remain_temp=[temp;NaN*(ones(MaxRep_temp*length(unique_data)-length(temp),1))];
        data_reshape(:,[i+1:end])=reshape(data_remain_temp,length(unique_data),MaxRep_temp);
        
        
    else
        RepData(:,i)=curr_data;
        
    end
    
end

NumInRep=sum(~isnan(RepData),1);

%Return the temporal location
temp_loc=NaN*ones(1,length(data));

for j=1:length(unique_data)
    curr_data=unique_data(j);
    [loc_curr,c,~]=find(RepData==curr_data);
    temp_loc(data==curr_data)=loc_curr';
end



end