function [reISI,reISISeq] = ISIInInterval(interval,time,bin);



MaxSampleBin = ceil(max(interval(:,2) -interval(:,1)));
reISI = NaN* ones(size(interval,1),MaxSampleBin);
MaxTimeBin = ceil(MaxSampleBin/bin);
reISISeq = NaN*ones(size(interval,1),MaxTimeBin);
time = time(:,1:end-2);


for i = 1:size(interval,1)
    val = interval(i,:);
    sel = time(i,:) < val(2) & time(i,:) >= val(1);
    spk_val = time(i,:);

    SampleBin = ceil(diff(val));
    timeBin = floor(SampleBin/bin);
    if sum(sel) >0
        
        isi_pre = diff(spk_val);
        %Remove those duplicate spikes(ISI < 1ms)
        dul = find(isi_pre< 1);
        if ~isempty(dul)
            spk_val(dul) =[];%Remove
            sel(dul) = [];
        end
        
        
        isi = 1./diff(spk_val)*1000;
        isi = [isi,NaN];
        fillzero = find(isnan(isi) == 1,1,'first');
        isi(fillzero) = 0;
        
       
        
        for j = 1:SampleBin
            interval_j = val(1) + j -1: val(1) + j;
            spk_time_j = spk_val< interval_j(2) & spk_val >= interval_j(1);
            loc = find(spk_time_j == 1,1,'first');
            
            
            if isempty(loc)
                %If no spike
                if j == 1
                    %
                    spkloc = find(sel == 1,1,'first')-1;%fill in with the last spike isi;
                    if spkloc == 0
                        preISI = 0;
                    else
                    
                        preISI = isi(spkloc);
                    
                    end
                    if isempty(preISI)
                        reISI(i,j) = NaN; 
                    else
                        reISI(i,j) = preISI; 
                    end
                    
                    
                    
                else
                    reISI(i,j) = reISI(i,j-1);
                end
            else
                %if there is spike, update 1/isi
                %
               
                reISI(i,j) = isi(loc);
                
                
             
            end
            
        end
        
        reISISeq_tmp = reshape(reISI(i,1:floor(SampleBin/bin)*bin),bin,floor(SampleBin/bin));
        reISISeq(i,1:timeBin) = mean(reISISeq_tmp,1);
        
        %Plot raster with 1/ISI
        %{
        StartLoc = 1;
        LineWidth = 2;
        Barlength = 1;
        spk_time_sel =spk_val(sel);
        close all;
        figure(1)
        
        hold on
        valbin =[val(1):MaxSampleBin+val(1)-1]+0.5;
        plot(valbin,reISI(i,:),'-k');
        SpikeRaterPlot(spk_time_sel,StartLoc,Barlength,'r',LineWidth);
        ylim([0,100]);
        %}

    end
    
    
    
    
end


end