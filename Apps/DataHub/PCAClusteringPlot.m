function PCAClusteringPlot(WaveformData,WaveformName,SpikeChannel,ShowFigureFlag,OutputFlag)
WaveformName=WaveformName(2:end-3);
%FileName=WaveformData.RawFileName;
%Order=WaveformData.Order;

%Reload Waveform Data
WaveformPath=WaveformData.Path;
DataType=WaveformPath(end-3:end);
if contains(DataType,'mat')
    WaveformData=WaveformPlexon(WaveformPath);
else
    %load from Blip original file
     WaveformData=LoadWaveform(WaveformPath);
    
end
PrepareFigure(1,'w',[100,100, 400,400],'Name','Principal Component Scatter Plot with Colored Clusters');   

%Plot the spikes in PCA coordinates; YXF 11082019
          SpikeType_All=length(unique(WaveformData.Type))-1;
            if SpikeChannel>3 & SpikeChannel<=3+SpikeType_All+1
                    SpikeType=SpikeChannel-3;
                elseif SpikeChannel<=3 
                    SpikeType=min(SpikeChannel,3);
                    
                elseif SpikeChannel==3+size(WaveformData.Wform,1)+1
                     SpikeType=1:min(SpikeChannel,3);
                elseif SpikeChannel==3+size(WaveformData.Wform,1)+2
                    SpikeType=1:SpikeType_All;
            end
            
            
            for i=1:length(SpikeType)
             
              subplot(ceil(length(SpikeType)/ceil(sqrt(length(SpikeType)))),ceil(sqrt(length(SpikeType))),i);

             [pc, zscores, pcvars] =pca(WaveformData.Wform(WaveformData.Type==SpikeType(i)-1,:));
             cumsum(pcvars./sum(pcvars) * 100);
             cla;legend('off');
             gscatter(zscores(:,1),zscores(:,2),i*ones(size(zscores(:,1),1),1),'rgbcym','.******');
             xlabel('PC1');ylabel('PC2');
           %  title(WaveformName);
             
             axis equal
           % plot(T0,WaveformData.Wform(WaveformData.Type==i-1,:)','color',[0.3+(i-1)/2,0.8-(i-1)/2,0.1])
           % xlim([-5000,5000]);
           % ylim([-5000,5000]);
            set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3)
            
            end
           
            


end