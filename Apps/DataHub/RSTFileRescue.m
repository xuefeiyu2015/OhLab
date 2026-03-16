function File=RSTFileRescue(RawFiles);

for f=1:length(RawFiles) %Sequentially go through and process each file.
          fname=RawFiles(f).name;    
          cd(RawFiles(f).folder);
          FID= fopen(fname);       
          % determining the binary file length
          fseek(FID,0,1);%0:Indicate reading offset;1: indicate start from end; Fseek sets the current location to the end of the file
          P_end=ftell(FID);%ftell returns the current reading position
          fseek(FID,0,-1); %Fseek sets current location back to the beginning.
          Trlnum=1; %Initiation of trial number
          clear RAW
          HdrLen=fread(FID,1,'int'); %read one interger from the binary file.The interger turned out to be HdrLen;
          %It seems that in the original binary file, the first integer
          RAW.FileName=fname;
          %indicates the following char number;
          RAW.Header=char(fread(FID,HdrLen,'char')'); % Next, read in HdrLen numbers of characters, which are the header.char returns 7-bit ASCII characters;After reading, the pointer will be moved to the end;
          ctr0=0; % Not sure.
              while ftell(FID)<P_end %while not at end of file
                  
                  RAW.Data{Trlnum}=fread(FID,1,'int');
                  
                  %{
                    RAW.TrlID(Trlnum,1)=fread(FID,1,'int'); %The the next interger is the Trial Index      
                    % The time output (in microSec) is converted to ms by 0.001* and made relative to the begining of the trial
                    RAW.TrlTime(Trlnum,:)=0.001*fread(FID,2,'int64')';%Then next, read in two numbers of  64 bit interger %Trial start tiem and end time(according to the event code?).
                    numEventCodes=fread(FID,1,'int'); %Number of event codes
                    if numEventCodes<0
                         RAW.numEventCodes(Trlnum,1)=NaN; %Number of event codes
                    RAW.EventCode{Trlnum,1}=NaN;
                    RAW.EventTime{Trlnum,1}=NaN;
                        
                    else
                    
                    RAW.numEventCodes(Trlnum,1)=numEventCodes; %Number of event codes
                    RAW.EventCode{Trlnum,1}=fread(FID,numEventCodes,'int')';
                    RAW.EventTime{Trlnum,1}=0.001*fread(FID,numEventCodes,'int64')';
        
                    % user added events , time here is also converted to ms 6/18
                    NumOfUserData=fread(FID,1,'int'); %What is Number of User Data?
                    WieldData=1;
                    %{
                    while NumOfUserData<0
                        RAW.NumOfUserData(Trlnum,1)=NumOfUserData;
                        NumOfUserData=fread(FID,1,'int');
                        DataHold(WieldData)=NumOfUserData;
                         WieldData= WieldData+1;
                    
                    end
                    %}
                    if NumOfUserData>0
                    RAW.NumOfUserData(Trlnum,1)=NumOfUserData;
                    
                         for ctr=ctr0+1:ctr0+NumOfUserData %What is ctr and ctr0?
                              SizeOfName=fread(FID,1,'int');
                              RAW.UserData{ctr,3}=char([fread(FID,SizeOfName,'char')]'); %ID,Switched,OnOff,X,Y
                              %ID is the object identity(e.g.fractal set number; Fixation point identity); Switched is whether
                              %it's point to the next process or else;
                              %OnOff indicated whether there is visual
                              %stimulus shown;
                              %X,Y indicates the location of the visual
                              %stimulus(location of the fixation point or the location of the target);
                              %Each userdata indicates different task procedural;
                              NumOfValues=abs(fread(FID,1,'int'));
                              RAW.UserData{ctr,1}=[Trlnum;fread(FID,NumOfValues,'double')]';%The values related to ID Switched, OnOff,X,Y
                              RAW.UserData{ctr,2}=[0.001*fread(FID,1,'int64')]; %Trial start time 
                         end
                    
        
                           % this line was added to deal with NumOfUserData=0, 9/7/12
                         if NumOfUserData>0, ctr0=ctr; end
       
                           % Spike Times
                            numChans=fread(FID,1,'int');
                            if numChans<=0
                                RAW.SpkTime{Trlnum,1}=[];
                            end
                            
                          
                         for chan=1:numChans
                             NChanID(chan)=fread(FID,1,'int');
                             numNns=fread(FID,1,'int');
                             for Nn=1:numNns
                              NnID(Nn)=fread(FID,1,'int');
                              numSpikesOfNn(Nn)=fread(FID,1,'int');
                              RAW.SpkTime{Trlnum,NChanID(chan)*3+NnID(Nn)+1}=[0.001*fread(FID,numSpikesOfNn(Nn),'int64')]';
                              end
                         end
        
                        
                           % Analoge Input
                           RAW.AnalogSRate=fread(FID,1,'int');%The sampling rate of the analog signal
                           RAW.AnalogTime(Trlnum,:)=[0.001*fread(FID,2,'int64')]';%The trial start time and the trial end time
            
                           numBehaviorChannels=fread(FID,1,'int');%Number of behavior channel saved
                           if datenum(fname(1:6),'mmddyy')<=datenum('022112','mmddyy') % changes in saving data after a particular date
                              numDataPointsOfBehaviorAI=fread(FID,1,'int');
                           end
                           for chan=1:numBehaviorChannels
                               AChanID(chan)=fread(FID,1,'int');%First two are Eye channels: 1:x,Horizontal; 2:y,Vertical;
                                if datenum(fname(1:6),'mmddyy')>datenum('022112','mmddyy')
                                   numDataPointsOfBehaviorAI=fread(FID,1,'int');
                                end
                                if numDataPointsOfBehaviorAI>0
                                   RAW.AnalogAmp{Trlnum,AChanID(chan)+1}=[fread(FID,numDataPointsOfBehaviorAI,'int')]';
                                else
                                   RAW.AnalogAmp{Trlnum,AChanID(chan)+1}=NaN;
                                end 
                           end
                        
                 % hold on; plot(RAW.AIAmp{Trlnum,1}(1,:),RAW.AIAmp{Trlnum,1}(2,:));
                    end
                    end
                  %}
              Trlnum=Trlnum+1;
              end


end