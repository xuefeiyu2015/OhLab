function ave = CalculateEveryNTrial(data,N,Step,tolN,med);
%Function to calculate the moving everage of data for N trials, smoothing
%by Step
% data: the original data sequence before smoothing; row: each trial, col:
% data for each time step
% N: Number of data point to do the average;
% Step: the moving step for the smoothing 
%TolN: Tolerance Number of trials to calculate the Mean in the end
% Output: ave: the resulting moving average for the dataset after the
% average;ave.data = data seq after smoothing; ave.time = resulting bin seq
col = size(data,2) -tolN + 1;
total_bin = floor(col/Step);
data_o = zeros(size(data,1),total_bin);
stp_idx = 0;
for i = 1:Step: col
    stp_idx = stp_idx + 1;
    if i+N-1 > col
        sel = data(:,i:end);
    else
        sel = data(:,i:i+N-1);
     
        
    end
    if med == 1
        data_o(:, stp_idx) = nanmedian(sel,2);
    else
    data_o(:, stp_idx) = nanmean(sel,2);
    end
    
    
    
    
end
ave.data = data_o;
ave.time = 1:Step:col;


end