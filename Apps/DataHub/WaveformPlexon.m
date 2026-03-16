function wv=WaveformPlexon(path)
%Function to load Waveform
%Load from Plexon
    load(path);
    SpikeWaveForm_Extra=wv_task.Waveform;
    SR=wv_task.SpikeWaveFormSamplingRate;
    
    %Packing the waveform sample to accomondate the old styles using
          %the median waveform
          WaveFormSampleNum=1000;
          NumOfFake=2;
          SpikeWaveFormSample_Wform=cell2mat(cellfun(@(x) x(1:min(size(x,1),WaveFormSampleNum),:),SpikeWaveForm_Extra,'uniform',0)');
          
         
          SpikeWaveFormSample_Num=cell2mat(cellfun(@(x) min(size(x,1),WaveFormSampleNum),SpikeWaveForm_Extra,'uniform',0)');
         SpikeWaveFormSample_Num_Array=1:sum(SpikeWaveFormSample_Num)+NumOfFake;
       %  SpikeWaveFormSample_Num_Array=1:sum(SpikeWaveFormSample_Num)+NumOfFake;
       
       %{
          SpikeWaveFormSample_Type=-1*ones(sum(SpikeWaveFormSample_Num)+NumOfFake,1); %Add some fake data about the unsorted spikes
          Num_tmp=cumsum(SpikeWaveFormSample_Num);
          Num_tmp=[0;Num_tmp];
         
          for i=1:max(length(SpikeWaveFormSample_Num)-1,1)
              %SpikeWaveFormSample_Type((i-1)*SpikeWaveFormSample_Num(max(i-1,1))+1:i*SpikeWaveFormSample_Num(i))=-1+i;
              SpikeWaveFormSample_Type(Num_tmp(i)+1:Num_tmp(i+1))=-1+i;
    
          end
          %}
          SpikeTypeNum=size(SpikeWaveForm_Extra,2);
           SpikeWaveFormSample_Type=[];
           for i=1:SpikeTypeNum
               SpikeWaveFormSample_Type=[SpikeWaveFormSample_Type;(-1+i)*ones(SpikeWaveFormSample_Num(i),1)];
               
           end
          SpikeWaveFormSample_Type=[SpikeWaveFormSample_Type;-1*ones(NumOfFake,1)];
          
          %Add some fake data about the unsorted spikes
         SpikeWaveFormSample_Wform=[SpikeWaveFormSample_Wform;NaN*ones(NumOfFake,size(SpikeWaveFormSample_Wform,2))];
         SampleSR=SR*ones(sum(SpikeWaveFormSample_Num)+NumOfFake,1);
         
         SpikeWaveFormSample.Wform=SpikeWaveFormSample_Wform;
         SpikeWaveFormSample.SR=SampleSR;
         SpikeWaveFormSample.Num=SpikeWaveFormSample_Num_Array;
         SpikeWaveFormSample.Type=SpikeWaveFormSample_Type;
         SpikeWaveFormSample.MedWave=cell2mat(cellfun(@(x) median(x,1,'omitnan' ), SpikeWaveForm_Extra,'uniform',0)');
         
         %}
wv=SpikeWaveFormSample;
end