function data=InsertNullTrials(data_ori);

%Function to insert null trials for whole data analysis
%by xuefei yu 05262021
Rownum=max(data_ori(:,1))-min(data_ori(:,1))+1;

data=NaN*ones(Rownum,size(data_ori,2));
data(:,1)=min(data_ori(:,1)):max(data_ori(:,1));
data(ismember(data(:,1),data_ori(:,1)),2:end)=data_ori(:,2:end);

end