function SpikeRaterPlot(SpikeTime,StartLoc,Barlength,Color,LineWidth)
%Plot spike raster here;YXF 10162019
Xplot=[];
Yplot=[];
  for i=1:size(SpikeTime,1)
      
      Xplot=[Xplot,repmat(SpikeTime(i,:),2,1)];
      Yplot=[Yplot,repmat([StartLoc+i-1;StartLoc+i-1+Barlength],1,size(SpikeTime(i,:),2))];
  
  end
  
  
  
plot(Xplot, Yplot,'Color',Color,'LineWidth',LineWidth);
ylim([Barlength,StartLoc-1+size(SpikeTime,1)+Barlength*2]);
%yticks ([ceil(Barlength):floor(size(SpikeTime,1)/2):floor(i+Barlength*2)])

end

