% load spike wave form assumes you have *.wv file saved in  \Neural Raw Data\ and will save in folder \Neural Mat Data
% to load a file modify line 8: WavFiles=dir([fpath,'043014SPK*.wv']); with the file name and run use wildcard to load and save multiple files
% Ali Ghazizadeh 9/14/2016
%Modefied by YXF 11012019
function WaveformData=LoadWaveform(path)
     WavFiles=dir(path);

% Made from DAS_cl, DAS_View.cpp function SaveWaves
for f=1:length(WavFiles)
    fname=WavFiles(f).name;
    cd(WavFiles(f).folder);
 
    [FID, MESSAGE] = fopen(fname); 
    % determining the binary file length
    fseek(FID,0,1);P_end=ftell(FID);fseek(FID,0,-1); Wvnum=1;clear WAV;
    
    while ftell(FID)<P_end
        WAV.Num(Wvnum,1)= fread(FID,1,'int'); % wav number to be saved
        WAV.SR(Wvnum,1)= fread(FID,1,'int'); % sampling rate Hz
        WAV.WTime(Wvnum,1)= 0.001*fread(FID,1,'int64'); % wave timestamp,  time converterd to ms from microsec
        WAV.WTime2Threshold(Wvnum,1)= fread(FID,1,'double');% The time of the first signal
        WAV.DataPointNum(Wvnum,1)= fread(FID,1,'int');
        WAV.Type(Wvnum,1)=fread(FID,1,'int'); % The ECode: not selected(<0), neuron 1 (0), neuron 2(1), neuron3(2).
        WAV.ChanId(Wvnum,1)=fread(FID,1,'int'); % recording channel ID
        WAV.Dummy(Wvnum,1)=fread(FID,1,'int'); % FutureUseInt
        WAV.Wform(Wvnum,:)=fread(FID,WAV.DataPointNum(Wvnum,1),'int'); % FutureUseInt
        Wvnum=Wvnum+1;
    end
    
    
     fclose(FID); % closing the file after saving 5/15/2013
     WaveformData=WAV;
end


end




% to view the first neuron spikes use
%T0=[0: unique(WAV.SR(WAV.Type==0,1)):(unique(WAV.DataPointNum(WAV.Type==0,1))-1)*unique(WAV.SR(WAV.Type==0,1))];
%plot(T0,WAV.Wform(WAV.Type==0,:)')