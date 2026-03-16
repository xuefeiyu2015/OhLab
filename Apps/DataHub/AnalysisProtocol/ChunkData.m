function ChunkedData=ChunkData(TrialIndex,Data,OnsetBin,OffsetBin)
%This function is to chunk data into small parts according to the start bin
%/time and the stop bin/time
%It is initially designed for the eyemovent analysis in bonus trials in the foraging task.
%Xuefei Yu, 02182020
ChunkedData_temp=ones(numel(OnsetBin),size(Data,2))*NaN;
maxLength=1;
maxNum=1;
for i=1:length(TrialIndex)
    CurrData=Data(i,:);
    
    CheckIndex=isnan(OnsetBin(i,:));
    CurrDataChunk=repmat(CurrData,sum(~CheckIndex),1);
    ChunkInterval=[OnsetBin(i,~CheckIndex)',OffsetBin(i,~CheckIndex)']; 
    CurrentChunk=SelectMatrixData(CurrDataChunk,ChunkInterval(:,1),ChunkInterval(:,2));
    CurrentChunk=[repmat(TrialIndex(i),sum(~CheckIndex),1),CurrentChunk];
    ChunkedData_temp(maxNum:maxNum+sum(~CheckIndex)-1,1:size(CurrentChunk,2))=CurrentChunk;
    maxNum=maxNum+sum(~CheckIndex);
    if maxLength<size(CurrentChunk,2)
        maxLength=size(CurrentChunk,2);
    end
end
ChunkedData=ChunkedData_temp(1:maxNum-1,1:maxLength);
end