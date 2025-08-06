function Batch_SequenceLearning(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);

BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
%{
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
%}
PlotExample = 1;
%Load file
FilesName_Each=Data.FileName;
ChannelID=Data.ChannelNumber;
ChannelFull=[];
index = 1;

seq_cr =  containers.Map;
seq_tt =  containers.Map;
seq_r = containers.Map;

seq_cri =  containers.Map;

seq_cr_each =  containers.Map;

seq_correctrate_bin =  containers.Map;

SeqCount = containers.Map;

seq_ctime =  containers.Map;

seq_ave_ctime = containers.Map;


seq_TimeInterval =  containers.Map;
seq_DateLast =  containers.Map;

seq_TrainingDays =  containers.Map;

PreGapTraining_ct=containers.Map;

Ave_Dist_All = containers.Map;
Distract_Ave_Dist_All=containers.Map;

Ave_Angle_All = containers.Map;
Distract_Ave_Angle_All=containers.Map;

r_seq_All = containers.Map;
p_seq_All=containers.Map;

seq_cri_post=  containers.Map;

%Monkey_save =[];
for i=1:length(FilesName)
    for j = 1:length(FilesName{i})
        clear OutputData;
        ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];
  
        load(FilesName{i}{j});
      %  disp(FilesName{i}{j})
     
     for n = 1:length(OutputData.SequenceLearning)
          Task(index)=OutputData.SequenceLearning(n).TaskCode;
          TaskName{index}=OutputData.SequenceLearning(n).Task;
          FilesNameAll{index}=FilesName_Each{i};

          Monkey = FilesNameAll{index}(1);

          if Monkey == 'R'
              MonkeyIndex = 1;
          else
              MonkeyIndex = 2;
          end
        %  Monkey_save(index) = MonkeyIndex;

          Date_curr = RecordDate{i};
          Data =OutputData.SequenceLearning(n).DataStamp;

          target_date = datetime('04/26/2022', 'InputFormat', 'MM/dd/yyyy');
          curr_date = Date_curr;

          previous = curr_date <target_date;

          %% Load the parameters

          %% Hyperset
          Seq = Data('UniqueSeq');

          if MonkeyIndex==1
              Seq=Seq+2000;
          
          end

          % Correct rate
          CorrectRate= Data('CorrectRate');

            TgtAcTime = Data('TgtAcTime');
            TotalCompleteTrial = Data('TotalCompleteTrial');

            CorrectTrial = Data('CorrectTrial');

             NTrialToCri= Data('NTrialToCri');

             TrialOver60 = Data('TrialOver60');

             r_seq = Data('r_seq');

            p_seq = Data('p_seq');

            %Time 

            TgtATimeEach_Mean = Data('TgtATimeEach_Mean');

            HypersetCTSeq = Data('HypersetCTSeq');
            HypersetCT_ave = Data('HypersetCT_ave');

            Hyperset_Ave_AT = Data('Hyperset_Ave_AT');
             Hyperset_Seq_Ave = Data('Hyperset_Seq_Ave');


            CorrectRate_seq = Data('LR_Seq');%,'LR_Bin'

            %Average distances
            Distract_Ave_Dist = Data('Dist_Seq_Dis_Mean');

            Ave_Dist = Data('Dist_Seq_Mean');

            %Average angle

             Distract_Ave_Angle = Data('Angle_Seq_Dis_Mean');

            Ave_Angle = Data('Angle_Seq_Mean');


            


            %% Sets in a hyperset
            TrialOver60_j = Data('TrialOver60_j');

            NTrialToCri_j = Data('NTrialToCri_j');

           

            Hyperset_Ave_AT_j = Data('Hyperset_Ave_AT_j');

            TgtATimeEach_Mean = Data('TgtATimeEach_Mean');
            TgtATimeEach_Sem = Data('TgtATimeEach_Sem');

          %  r_seq = Data('r_seq');
          %  p_seq = Data('p_seq');


         %% Store completed sequence 

         % Correct Rate
         sel_store = CorrectTrial > 10;%Only store those which has been completed for 10 trials

         if sum(sel_store)>0

         seq_store = Seq(sel_store);
          CorrectRate_store = CorrectRate(sel_store);

          NTrialToCri_store = NTrialToCri(sel_store);

          

         % HypersetCTSeq_store = HypersetCTSeq(sel_store);
          HypersetCTSeq_store = Hyperset_Ave_AT(sel_store);
          
          
          CorrectRate_seq_store =CorrectRate_seq(sel_store);
          HypersetCT_ave_store =Hyperset_Seq_Ave(sel_store);

           AveDistance_seq_store =Ave_Dist(sel_store);
           Distract_AveDistance_seq_store =Distract_Ave_Dist(sel_store);


           AveAngle_seq_store =Ave_Angle(sel_store);
           Distract_AveAngle_seq_store =Distract_Ave_Angle(sel_store);

           r_seq_store =r_seq(sel_store);
           p_seq_store =p_seq(sel_store);


           %r_seq_All



 %          Ave_Angle_All = containers.Map;
%Distract_Ave_Angle_All=containers.Map;


         for m = 1:length(seq_store)
             

            

               %Sequence Training Date
             if isKey(seq_DateLast,string(seq_store(m)))
                 seq_DateLast(string(seq_store(m)))=[seq_DateLast(string(seq_store(m))),{Date_curr}];
             else
                 seq_DateLast(string(seq_store(m)))= {Date_curr};
             end

             %Sequence Training Date Interval
             if isKey(seq_TimeInterval,string(seq_store(m)))
                 lastdates =  seq_DateLast(string(seq_store(m)));
                 interval_dates = days(datetime(Date_curr) - datetime(lastdates{end-1}));
                 


                 seq_TimeInterval(string(seq_store(m)))=[seq_TimeInterval(string(seq_store(m))),days(datetime(Date_curr) - datetime(lastdates{end-1}))];
             else
                 seq_TimeInterval(string(seq_store(m)))= 0;
                 interval_dates = 0;
             end
%
            
                 
%}
%GapAdd=0;%Remove this 

            %Correct Rate

            if isKey(seq_cr,string(seq_store(m)))
                
                
                seq_cr(string(seq_store(m))) = [seq_cr(string(seq_store(m))),CorrectRate_store(m)];
               
             
            else

                seq_cr(string(seq_store(m))) = CorrectRate_store(m);
            end

            %Number of trials to criterium

             if isKey(seq_cri,string(seq_store(m)))
                 
                    seq_cri(string(seq_store(m))) = [seq_cri(string(seq_store(m))),NTrialToCri_store(m)];
                
             
            else

                seq_cri(string(seq_store(m))) = NTrialToCri_store(m);
             end


             %Correct rate sequence

             if isKey(seq_correctrate_bin,string(seq_store(m)))
               
                    seq_correctrate_bin(string(seq_store(m))) = [seq_correctrate_bin(string(seq_store(m))),CorrectRate_seq_store(m)];
                
             
            else

                seq_correctrate_bin(string(seq_store(m))) = CorrectRate_seq_store(m);
             end

             %Average distance 
              if isKey(Ave_Dist_All,string(seq_store(m)))
                
                    Ave_Dist_All(string(seq_store(m))) = [Ave_Dist_All(string(seq_store(m))),AveDistance_seq_store(m)];
               
             
            else

                Ave_Dist_All(string(seq_store(m))) = AveDistance_seq_store(m);
              end

              %Distract_Ave_Dist_All
              if isKey(Distract_Ave_Dist_All,string(seq_store(m)))
                
                     Distract_Ave_Dist_All(string(seq_store(m))) = [ Distract_Ave_Dist_All(string(seq_store(m))), Distract_AveDistance_seq_store(m)];
                
             
            else

                 Distract_Ave_Dist_All(string(seq_store(m))) =  Distract_AveDistance_seq_store(m);
              end

              %Average angle
              if isKey(Ave_Angle_All,string(seq_store(m)))
                
                    Ave_Angle_All(string(seq_store(m))) = [Ave_Angle_All(string(seq_store(m))),AveAngle_seq_store(m)];
                
             
            else

                Ave_Angle_All(string(seq_store(m))) = AveAngle_seq_store(m);
              end

              %Distract_Ave_Dist_All
              if isKey(Distract_Ave_Angle_All,string(seq_store(m)))
               
                     Distract_Ave_Angle_All(string(seq_store(m))) = [ Distract_Ave_Angle_All(string(seq_store(m))), Distract_AveAngle_seq_store(m)];
                
             
            else

                 Distract_Ave_Angle_All(string(seq_store(m))) =  Distract_AveAngle_seq_store(m);
              end

               %r_seq_all
              if isKey(r_seq_All,string(seq_store(m)))
               
                     r_seq_All(string(seq_store(m))) = [ r_seq_All(string(seq_store(m))), r_seq_store(m)];
               
             
            else

                 r_seq_All(string(seq_store(m))) = r_seq_store(m);

              end

%p_seq_all
               if isKey(p_seq_All,string(seq_store(m)))
                
                     p_seq_All(string(seq_store(m))) = [ p_seq_All(string(seq_store(m))), p_seq_store(m)];
               
             
            else

                 p_seq_All(string(seq_store(m))) =  p_seq_store(m);
              end



             

             

             % Sequence number count
            if isKey(SeqCount,string(seq_store(m)))
                
               SeqCount(string(seq_store(m))) = SeqCount(string(seq_store(m))) +1;
             
            else

               SeqCount(string(seq_store(m))) = 1;
            end

            %Sequence completion time

            if isKey(seq_ctime,string(seq_store(m)))

                            
                    seq_ctime(string(seq_store(m))) = [seq_ctime(string(seq_store(m))),{HypersetCTSeq_store{m}*5}];
               
             
            else

               seq_ctime(string(seq_store(m))) = {HypersetCTSeq_store{m}*5};
            end

            %Average sequence completion time
             if isKey(seq_ave_ctime,string(seq_store(m)))
                 
                    seq_ave_ctime(string(seq_store(m))) = [seq_ave_ctime(string(seq_store(m))),HypersetCT_ave_store(m)*5];
                 
             else
                 seq_ave_ctime(string(seq_store(m))) = HypersetCT_ave_store(m)*5;
             end

           

          if previous==0

              % else %After the gap

                     if isKey(seq_cri_post,string(seq_store(m)))
                
                
                        seq_cri_post(string(seq_store(m))) = [seq_cri_post(string(seq_store(m))),NTrialToCri_store(m)];
                
             
                    else

                        seq_cri_post(string(seq_store(m))) = NTrialToCri_store(m);
                     end
             



  %  seq_cri_post = 

         end %End of if previous ==1


         end %End of if length(seq_store)

           
         
         end %End of if isempty


%{

            'UniqueSeq','CorrectRate','TgtAcTime','TotalCompleteTrial','CorrectTrial','NTrialToCri','TrialOver60',...
    'LR_Seq','LR_Bin',...
    'TrialOver60_j','NTrialToCri_j',...
    'HypersetCTSeq','HypersetCT_ave','Hyperset_Ave_AT','Hyperset_Seq_Ave',...
    'Hyperset_Ave_AT_j',...
    'r_seq','p_seq','TgtATimeEach_Mean','TgtATimeEach_Sem','Dist_Seq_Mean','Dist_Seq_Sem','Angle_Seq_Mean','Angle_Seq_Sem',...
    'Dist_Seq_Dis_Mean','Dist_Seq_Dis_Sem','Angle_Seq_Dis_Mean','Angle_Seq_Dis_Sem'
    
%}

     end

    

    end
end

seq_all = keys(seq_correctrate_bin);

%Correct rate
all_data_correct_rate_sequence = seq_correctrate_bin.values;

%Completion time seq

all_data_ctime_sequence = seq_ctime.values;

%Average completion time
all_data_ctime_ave = seq_ave_ctime.values;
all_data_ctime_ave=cellfun(@(x) x',all_data_ctime_ave,'UniformOutput',0);
all_data_ctime_ave_mat = fillInNaN(all_data_ctime_ave);

%Before gap training days 
%training_days = seq_TrainingDays.values;
%training_gap_seqs = keys(seq_TrainingDays);

%all_data_ctime_beforegap = cellfun(@(x) x',all_data_ctime_ave,'UniformOutput',0);

%{
PreGapTraining_ct_val = PreGapTraining_ct.values;
PreGapTraining_ct_val  = cellfun(@(x) x',PreGapTraining_ct_val,'UniformOutput',0);
PreGapTraining_ct_mat= fillInNaN(PreGapTraining_ct_val);
%}


%dAve_Dist_All

Ave_Dist_All_val =Ave_Dist_All.values;
Ave_Dist_All_val_mat  = cellfun(@(x) x',Ave_Dist_All_val,'UniformOutput',0);
Ave_Dist_All_mat = fillInNaN(Ave_Dist_All_val_mat );

Dis_Ave_Dist_All_val = Distract_Ave_Dist_All.values;
Dis_Ave_Dist_All_val_mat =cellfun(@(x) x',Dis_Ave_Dist_All_val,'UniformOutput',0);
Dis_Ave_Dist_All_mat = fillInNaN(Dis_Ave_Dist_All_val_mat );



%dAve_Angle_All

Ave_Angle_All_val =Ave_Angle_All.values;
Ave_Angle_All_val_mat  = cellfun(@(x) x',Ave_Angle_All_val,'UniformOutput',0);
Ave_Angle_All_mat = fillInNaN(Ave_Angle_All_val_mat );

Dis_Ave_Angle_All_val = Distract_Ave_Angle_All.values;
Dis_Ave_Angle_All_val_mat =cellfun(@(x) x',Dis_Ave_Angle_All_val,'UniformOutput',0);
Dis_Ave_Angle_All_mat = fillInNaN(Dis_Ave_Angle_All_val_mat );


%r_seq and p_seq

r_seq_val =r_seq_All.values;
r_seq_val_mat  = cellfun(@(x) x',r_seq_val,'UniformOutput',0);
r_seq_mat = fillInNaN(r_seq_val_mat );

p_seq_val = p_seq_All.values;
p_seq_val_mat  =cellfun(@(x) x',p_seq_val,'UniformOutput',0);
p_seq_mat = fillInNaN(p_seq_val_mat );




%% Initial training
Initial_learning_cr = cellfun(@(x) x{1}',all_data_correct_rate_sequence,'UniformOutput',0);

InitialLearningCurve= fillInNaN(Initial_learning_cr);


tmp = seq_cri.values;

NTrialCri_Initial =  cellfun(@(x) x(1),tmp );


SeqCount_Val = cell2mat(SeqCount.values);


Initial_learning_ct = cellfun(@(x) x{1}',all_data_ctime_sequence,'UniformOutput',0);
InitialLearning_CompleteTime= fillInNaN(Initial_learning_ct);

%InitialLearning_CT_Ave = cellfun(@(x) x(1),all_data_ctime_ave);
InitialLearning_CT_Ave = all_data_ctime_ave_mat(:,1);


%Ave_Dist_All_mean = mean(Ave_Dist_All_mat(:,1));

r_seq_initial = r_seq_mat(:,1);
p_seq_initial = p_seq_mat(:,1);


%Compare the NToCri between scene based sequence and fractal based sequence

seq_trans =cellfun(@(x) str2num(x),seq_all);

frac_sel = seq_trans <1000 | (seq_trans >=2000 & seq_trans <3000);
scene_sel = seq_trans >= 1000| seq_trans >=3000;

NTrialCri_Initial_frac = NTrialCri_Initial(frac_sel);
NTrialCri_Initial_scene = NTrialCri_Initial(scene_sel);

Mean_NCri_frac = mean(NTrialCri_Initial_frac);
Mean_NCri_scene = mean(NTrialCri_Initial_scene);

disp('Frac N = ');
disp(sum(frac_sel));
disp('Scene N = ');
disp(sum(scene_sel));




%Compare the Completion time between scene based sequence and fractal based sequence
%Valid = InitialLearning_CT_Ave<=2000;
InitialLearning_CompleteTime_frac =InitialLearning_CT_Ave(frac_sel);%&Valid);
InitialLearning_CompleteTime_scene =InitialLearning_CT_Ave(scene_sel);%&Valid);

%InitialLearning_CompleteTime_frac=InitialLearning_CompleteTime_frac(InitialLearning_CompleteTime_frac<1000);
%InitialLearning_CompleteTime_scene=InitialLearning_CompleteTime_scene(InitialLearning_CompleteTime_scene<1000);


Mean_CT_frac = mean(InitialLearning_CompleteTime_frac);
Mean_CT_scene = mean(InitialLearning_CompleteTime_scene);


%Compare between pregap learning and post gap learning
%seq_cri_post

intervals = seq_TimeInterval.values;

has_gap = cellfun(@(x) sum(x>=30),intervals);

Gap_seq = seq_trans(has_gap>=1);




NTrialCri_Initial_before_gap = NTrialCri_Initial(has_gap>=1);


NTrialCri_post =seq_cri_post.values;
NTrialCri_post_initial_all= cellfun(@(x) x(1),NTrialCri_post);



Seq_post = keys(seq_cri_post);
Seq_post_mat = cellfun(@(x) str2num(x), Seq_post);

Gap_seq_select_post = ismember(Seq_post_mat,Gap_seq );


NTrialCri_post_initial_gap = NTrialCri_post_initial_all(Gap_seq_select_post);

%Find out the index of the gap
index_gap_all = cellfun(@(x) max(find(x>=30,1,'first')-1,1),intervals,'UniformOutput',0);
index_gap = index_gap_all(has_gap>=1);

tmp_gap = tmp(has_gap>=1);
NTrialCri_last_before_gap = cellfun(@(x,y) x(y),tmp_gap,index_gap);

% Comparison of the completion time before and after the gap
all_data_ctime_ave_gap = all_data_ctime_ave(has_gap>=1);

ctime_last_before_gap = cellfun(@(x,y) x(y),all_data_ctime_ave_gap ,index_gap);

index_aftergap = cellfun(@(x) x+1,index_gap,'UniformOutput',0); 

ctime_last_after_gap = cellfun(@(x,y) x(y),all_data_ctime_ave_gap ,index_aftergap);


c_time_initial_gap = InitialLearning_CompleteTime(has_gap>=1);
disp('Gap comparison N=?')
disp(length(c_time_initial_gap));

%{
 selec_valid = c_time_initial_gap<4000;

ctime_last_before_gap=ctime_last_before_gap(selec_valid);
ctime_last_after_gap = ctime_last_after_gap(selec_valid);
c_time_initial_gap = c_time_initial_gap(selec_valid);

NTrialCri_last_before_gap=NTrialCri_last_before_gap(selec_valid);
NTrialCri_post_initial_gap =NTrialCri_post_initial_gap(selec_valid);
NTrialCri_Initial_before_gap =NTrialCri_Initial_before_gap(selec_valid);
%}





%% Trained over 5 days 

WellTrained5 = SeqCount_Val >= 5;%Over 5 days

welltrained_tmp = tmp(WellTrained5);

NTrialCri_WellTrained5 =  cellfun(@(x) x(5),welltrained_tmp);

NTrialCri_Initial_match5 = NTrialCri_Initial(WellTrained5);

learningcurve_trained5 = all_data_correct_rate_sequence(WellTrained5);
learning5_cr = cellfun(@(x) x{5}',learningcurve_trained5,'UniformOutput',0);

After5LearningCurve= fillInNaN(learning5_cr);


%maxedge = max([NTrialCri_Initial_match5,NTrialCri_WellTrained5]);

seq_trained5 = seq_all(WellTrained5);

ct_trained5 = all_data_ctime_sequence(WellTrained5);
learning5_ct = cellfun(@(x) x{5}',ct_trained5,'UniformOutput',0);

After5LearningCT= fillInNaN(learning5_ct);

%ct_ave5 = all_data_ctime_ave(WellTrained5);

%cr_ave_ave5 = cellfun(@(x) x(5),ct_ave5,'UniformOutput',0);
InitialLearning_CT_Ave5 = all_data_ctime_ave_mat(:,5);
InitialLearning_CT_Ave5 = InitialLearning_CT_Ave5(~isnan(InitialLearning_CT_Ave5));




WellTrainedDays = 10;
WellTrained20 = SeqCount_Val >= WellTrainedDays;%Over 15 days
seq_trained20 = seq_all(WellTrained20);
welltrained_tmp = tmp(WellTrained20);

NTrialCri_WellTrained20 =  cellfun(@(x) x(WellTrainedDays),welltrained_tmp);

NTrialCri_Initial_match20 = NTrialCri_Initial(WellTrained20);

learningcurve_trained20 = all_data_correct_rate_sequence(WellTrained20);
learning20_cr = cellfun(@(x) x{WellTrainedDays}',learningcurve_trained20,'UniformOutput',0);

After20LearningCurve= fillInNaN(learning20_cr);

ct_trained20 = all_data_ctime_sequence(WellTrained20);
learning20_ct = cellfun(@(x) x{WellTrainedDays}',ct_trained20,'UniformOutput',0);

After20LearningCT= fillInNaN(learning20_ct);

InitialLearning_CT_Ave20 = all_data_ctime_ave_mat(:,WellTrainedDays);
InitialLearning_CT_Ave20 = InitialLearning_CT_Ave20(~isnan(InitialLearning_CT_Ave20));


r_seq_20 = r_seq_mat(:,WellTrainedDays);
p_seq_20 = p_seq_mat(:,WellTrainedDays);

r_seq_20 = r_seq_20(~isnan(r_seq_20));
p_seq_20 = p_seq_20(~isnan(p_seq_20));


%Completion time with days 
MoreThanOne = SeqCount_Val >1;

try
Ave_ct_sve_pretrain = nanmean(all_data_ctime_ave_mat(MoreThanOne,1:20),1);
Std_ct_sve_pretrain = nanstd(all_data_ctime_ave_mat(MoreThanOne,1:20),[],1)./sqrt(sum(~isnan(all_data_ctime_ave_mat(MoreThanOne,1:20)),1));
catch
    Ave_ct_sve_pretrain = nanmean(all_data_ctime_ave_mat(MoreThanOne,1:15),1);
Std_ct_sve_pretrain = nanstd(all_data_ctime_ave_mat(MoreThanOne,1:15),[],1)./sqrt(sum(~isnan(all_data_ctime_ave_mat(MoreThanOne,1:15)),1));
end

%Separate for each monkey
Monkey_save = ones(size(seq_trans));
Monkey_save(seq_trans <2000)=2;

Ave_ct_sve_pretrain_Robin = nanmean(all_data_ctime_ave_mat(MoreThanOne & Monkey_save==1,1:13),1);
Std_ct_sve_pretrain_Robin = nanstd(all_data_ctime_ave_mat(MoreThanOne & Monkey_save==1,1:13),[],1)./sqrt(sum(~isnan(all_data_ctime_ave_mat(MoreThanOne& Monkey_save==1,1:13)),1));
try
Ave_ct_sve_pretrain_Adams = nanmean(all_data_ctime_ave_mat(MoreThanOne & Monkey_save==2,1:20),1);
Std_ct_sve_pretrain_Adams = nanstd(all_data_ctime_ave_mat(MoreThanOne & Monkey_save==2,1:20),[],1)./sqrt(sum(~isnan(all_data_ctime_ave_mat(MoreThanOne& Monkey_save==2,1:20)),1));
catch
    Ave_ct_sve_pretrain_Adams = NaN*Ave_ct_sve_pretrain_Robin;
    Std_ct_sve_pretrain_Adams = NaN*Std_ct_sve_pretrain_Robin;

end



try

Ave_Dist_All_mat_mean = nanmean(Ave_Dist_All_mat(MoreThanOne,1:20),1);
Dis_Ave_Dist_All_mat_mean = nanmean(Dis_Ave_Dist_All_mat(MoreThanOne,1:20),1);

Ave_Dist_All_mat_std = nanstd(Ave_Dist_All_mat(MoreThanOne,1:20),[],1)./sqrt(sum(~isnan(Ave_Dist_All_mat(MoreThanOne,1:20)),1));
Dis_Ave_Dist_All_mat_std = nanstd(Dis_Ave_Dist_All_mat(MoreThanOne,1:20),[],1)./sqrt(sum(~isnan(Dis_Ave_Dist_All_mat(MoreThanOne,1:20)),1));
catch
    Ave_Dist_All_mat_mean = nanmean(Ave_Dist_All_mat(MoreThanOne,1:15),1);
Dis_Ave_Dist_All_mat_mean = nanmean(Dis_Ave_Dist_All_mat(MoreThanOne,1:15),1);

Ave_Dist_All_mat_std = nanstd(Ave_Dist_All_mat(MoreThanOne,1:15),[],1)./sqrt(sum(~isnan(Ave_Dist_All_mat(MoreThanOne,1:15)),1));
Dis_Ave_Dist_All_mat_std = nanstd(Dis_Ave_Dist_All_mat(MoreThanOne,1:15),[],1)./sqrt(sum(~isnan(Dis_Ave_Dist_All_mat(MoreThanOne,1:15)),1));
end

try
Ave_Angle_All_mat_mean = nanmean(Ave_Angle_All_mat(MoreThanOne,1:20),1);
Dis_Ave_Angle_All_mat_mean = nanmean(Dis_Ave_Angle_All_mat(MoreThanOne,1:20),1);

Ave_Angle_All_mat_std = nanstd(Ave_Angle_All_mat(MoreThanOne,1:20),1)./sqrt(sum(~isnan(Ave_Angle_All_mat(MoreThanOne,1:20)),1));
Dis_Ave_Angle_All_mat_std = nanstd(Dis_Ave_Angle_All_mat(MoreThanOne,1:20),1)./sqrt(sum(~isnan(Dis_Ave_Angle_All_mat(MoreThanOne,1:20)),1));
catch
    Ave_Angle_All_mat_mean = nanmean(Ave_Angle_All_mat(MoreThanOne,1:15),1);
Dis_Ave_Angle_All_mat_mean = nanmean(Dis_Ave_Angle_All_mat(MoreThanOne,1:15),1);

Ave_Angle_All_mat_std = nanstd(Ave_Angle_All_mat(MoreThanOne,1:15),1)./sqrt(sum(isnan(Ave_Angle_All_mat(MoreThanOne,1:15)),1));
Dis_Ave_Angle_All_mat_std = nanstd(Dis_Ave_Angle_All_mat(MoreThanOne,1:15),1)./sqrt(sum(~isnan(Dis_Ave_Angle_All_mat(MoreThanOne,1:15)),1));
end
try
r_seq_mean = nanmean(r_seq_mat(:,1:20),1);
r_seq_sem = nanstd(r_seq_mat(:,1:20),[],1)./sqrt(sum(~isnan(r_seq_mat(:,1:20)),1));
catch
    r_seq_mean = nanmean(r_seq_mat(:,1:15),1);
r_seq_sem = nanstd(r_seq_mat(:,1:15),[],1)./sqrt(sum(~isnan(r_seq_mat(:,1:15)),1));
end


Example_LC = {'1007','1002','34','38'};
%Example_LC = {'1004','2'};
%Example_LC = {'1005'};



if ShowFigureFlag
 
FigureStartNum=100;
FigureIndex=1;

if PlotExample

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='Example_LearningCurve';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
Color=[Color;Color*0.5];
Color=[Color;Color*0.8];



subplot(3,1,1);
% Initial Learning Curve
for i = 1:length(Example_LC)
seq_curr_tmp = Example_LC{i};
sel = strcmp(seq_all, seq_curr_tmp);

ExampleLearningCurve = InitialLearningCurve(sel,:);


plot(1:length(ExampleLearningCurve),ExampleLearningCurve,'-b','color',Color(i,:),'LineWidth',3);
hold on

end
Example_LC_str = cellfun(@(x) sprintf('Seq: %s',x),Example_LC,'UniformOutput',0);
legend(Example_LC_str);
xlim([0,500]);
box off;

title('First day training');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

subplot(3,1,2);
% Training for 5 days
for i = 1:length(Example_LC)
seq_curr_tmp = Example_LC{i};
sel = strcmp(seq_trained5, seq_curr_tmp);

ExampleLearningCurve = After5LearningCurve(sel,:);


plot(1:length(ExampleLearningCurve),ExampleLearningCurve,'-b','color',Color(i,:),'LineWidth',3);
hold on

end

title('Trained for 5 days');

xlim([0,500]);
box off;
ylabel('Hyperset Completion Rate(%)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

subplot(3,1,3);
% Training for 20 days
for i = 1:length(Example_LC)
seq_curr_tmp = Example_LC{i};
sel = strcmp(seq_trained20, seq_curr_tmp);

ExampleLearningCurve = After20LearningCurve(sel,:);


plot(1:length(ExampleLearningCurve),ExampleLearningCurve,'-b','color',Color(i,:),'LineWidth',3);
hold on

end

title('Trained for 15 days');
xlabel('Trial Number');
xlim([0,500]);
box off;
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

end %End of if plot example




%% NTrialCri_Initial

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='NumberOfTrialsToCri';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 500,600],'Name',figtitlestr{FigureIndex});

all_num = [NTrialCri_Initial,NTrialCri_WellTrained5,NTrialCri_WellTrained20];
max_num = max(all_num);


edges = 1:10:max_num+10;

%Initial Learning
[N, edges] = histcounts(NTrialCri_Initial, edges,'Normalization', 'probability');

subplot(3,1,1);
bar(edges(1:end-1),N,'Facecolor','k','Edgecolor','k');
xlim([min(edges)-10,max(edges)]);
box off;

title('First day training');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');


%After 5
[N, edges] = histcounts(NTrialCri_WellTrained5, edges,'Normalization', 'probability');

subplot(3,1,2);
bar(edges(1:end-1),N,'Facecolor','k','Edgecolor','k');
xlim([min(edges)-10,max(edges)]);

title('Trained for 5 days');

box off;
ylabel('Proportion of Sequences');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');


%After 20
[N, edges] = histcounts(NTrialCri_WellTrained20, edges,'Normalization', 'probability');

subplot(3,1,3);
bar(edges(1:end-1),N,'Facecolor','k','Edgecolor','k');
xlim([min(edges)-10,max(edges)]);
xlabel('Number of Trials to Criterium');

box off;
title('Trained for 15 days');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');


if PlotExample
    %% Completion Time
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='CompletionTime';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 500,600],'Name',figtitlestr{FigureIndex});


% Initial Learning Completion Time
subplot(3,1,1);
for i = 1:length(Example_LC)
seq_curr_tmp = Example_LC{i};
sel = strcmp(seq_all, seq_curr_tmp);

ExampleLearningTime = InitialLearning_CompleteTime(sel,:);


plot(1:length(ExampleLearningTime),ExampleLearningTime,'-b','color',Color(i,:),'LineWidth',3);
hold on

end
Example_LC_str = cellfun(@(x) sprintf('Seq: %s',x),Example_LC,'UniformOutput',0);
legend(Example_LC_str);
%xlim([0,500]);
box off;

title('First day training');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

% After 5 days Completion Time
subplot(3,1,2);
for i = 1:length(Example_LC)
seq_curr_tmp = Example_LC{i};
sel = strcmp(seq_trained5, seq_curr_tmp);

ExampleLearningTime = After5LearningCT(sel,:);


plot(1:length(ExampleLearningTime),ExampleLearningTime,'-b','color',Color(i,:),'LineWidth',3);
hold on

end
%Example_LC_str = cellfun(@(x) sprintf('Seq: %s',x),Example_LC,'UniformOutput',0);
%legend(Example_LC_str);
%xlim([0,500]);
box off;
ylabel('Hyperset Completion Time(ms)');

title('Trained for 5 days');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

subplot(3,1,3);
% Training for 15 days
for i = 1:length(Example_LC)
seq_curr_tmp = Example_LC{i};
sel = strcmp(seq_trained20, seq_curr_tmp);

ExampleLearningCurve = After20LearningCT(sel,:);


plot(1:length(ExampleLearningCurve),ExampleLearningCurve,'-b','color',Color(i,:),'LineWidth',3);
hold on

end

title('Trained for 15 days');
xlabel('Trial Number');
%xlim([0,500]);
box off;
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');






%% Completion Time as days population
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PopulationCompletionTimeAsTraining';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,600],'Name',figtitlestr{FigureIndex});

all_num = [InitialLearning_CT_Ave;InitialLearning_CT_Ave5;InitialLearning_CT_Ave20];
max_num = max(all_num);


edges = 1:20:max_num+20;

%Initial Learning
[N, edges] = histcounts(InitialLearning_CT_Ave , edges,'Normalization', 'probability');

subplot(3,1,1);
bar(edges(1:end-1),N,'Facecolor','r','Edgecolor','r');
xlim([min(edges)-10,max(edges)]);
box off;

title('First day training');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');


%After 5
[N, edges] = histcounts(InitialLearning_CT_Ave5, edges,'Normalization', 'probability');

subplot(3,1,2);
bar(edges(1:end-1),N,'Facecolor','r','Edgecolor','r');
xlim([min(edges)-10,max(edges)]);

title('Trained for 5 days');

box off;
ylabel('Proportion of Sequences');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');


%After 20
[N, edges] = histcounts(InitialLearning_CT_Ave20, edges,'Normalization', 'probability');

subplot(3,1,3);
bar(edges(1:end-1),N,'Facecolor','r','Edgecolor','r');
xlim([min(edges)-10,max(edges)]);
xlabel('Hyperset Completion Time(ms)');

box off;
title('Trained for 15 days');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');


%errorbar(1:30,Ave_ct_sve_pretrain,Std_ct_sve_pretrain,'-k','LineWidth',3);

%% Completion Time as days example curve
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='ExampleCompletionTimeAsTraining';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,300],'Name',figtitlestr{FigureIndex});


set(fig{FigureIndex}, 'PaperUnits', 'inches');


% Completion time along with days 
for i = 1:length(Example_LC)
seq_curr_tmp = Example_LC{i};
sel = strcmp(seq_all, seq_curr_tmp);

ExampleLearningCurve = all_data_ctime_ave_mat(sel,1:20);



plot(1:length(ExampleLearningCurve),ExampleLearningCurve,'-b','color',Color(i,:),'LineWidth',3);
hold on

end
Example_LC_str = cellfun(@(x) sprintf('Seq: %s',x),Example_LC,'UniformOutput',0);
legend(Example_LC_str);

title('Comletion Time Along with Training');
xlabel('Days');
ylabel('Hyperset completion time(ms)');
%xlim([0,500]);
box off;
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

%for output
%{
set(fig{FigureIndex}, 'PaperUnits', 'inches');
formatFig(fig{FigureIndex}, [2 2], 'nature');

  FolderName='/Users/yux8/WorkRelated_YXF/Program_Matlab_Local/DataAnalysis/DataFromSystem/PackedData/SequenceFigures';
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
   
    OutputFigureName=strcat('Training','.pdf');
     saveas(fig{FigureIndex},OutputFigureName);
    disp('Figures have been exported to the neuron folder');

%}



end %End of if plot example

%% Completion Time as days population curve
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PopulationCompletionTimeAsTraining';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,300],'Name',figtitlestr{FigureIndex});
try
errorbar(1:20,Ave_ct_sve_pretrain,Std_ct_sve_pretrain,'-k','LineWidth',3);
catch
    errorbar(1:15,Ave_ct_sve_pretrain,Std_ct_sve_pretrain,'-k','LineWidth',3);
end

xlabel('Days');
ylabel('Hyperset completion time(ms)');
box off;
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

%% Completion Time as days population curve (separate monkey)
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PopulationCompletionTimeAsTrainingForEachMonkey';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,300],'Name',figtitlestr{FigureIndex});


errorbar(1:length(Ave_ct_sve_pretrain_Robin),Ave_ct_sve_pretrain_Robin,Std_ct_sve_pretrain_Robin,'-b','LineWidth',3);
hold on
errorbar(1:length(Ave_ct_sve_pretrain_Adams),Ave_ct_sve_pretrain_Adams,Std_ct_sve_pretrain_Adams,'-r','LineWidth',3);


xlabel('Days');
ylabel('Hyperset completion time(ms)');
legend({'Monkey R','Monkey A'});
box off;
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');




%Ave_ct_sve_pretrain_Robin

%% Distance and angle as days population curve
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='DistanceAngleAsTraining';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,300],'Name',figtitlestr{FigureIndex});
subplot(1,2,1);
plot(Ave_Dist_All_mat_mean,'-r','LineWidth',3);
hold on
errorbar(1:length(Ave_Dist_All_mat_mean),Ave_Dist_All_mat_mean,Ave_Dist_All_mat_std,'-r','LineWidth',3);

plot(Dis_Ave_Dist_All_mat_mean,'-b','LineWidth',3);
hold on
errorbar(1:length(Dis_Ave_Dist_All_mat_mean),Dis_Ave_Dist_All_mat_mean,Dis_Ave_Dist_All_mat_std,'-b','LineWidth',3);

box off;

xlabel('Training days');
ylabel('Average distance');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2);
plot(Ave_Angle_All_mat_mean,'-r','LineWidth',3);
hold on
errorbar(1:length(Ave_Angle_All_mat_mean),Ave_Angle_All_mat_mean,Ave_Angle_All_mat_std,'-r','LineWidth',3);


hold on
plot(Dis_Ave_Angle_All_mat_mean,'-b','LineWidth',3);
hold on
errorbar(1:length(Dis_Ave_Angle_All_mat_mean),Dis_Ave_Angle_All_mat_mean,Dis_Ave_Angle_All_mat_std,'-b','LineWidth',3);

box off;

xlabel('Training days');
ylabel('Average angle difference');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



%% r_seq histogram and change with days 
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='r_seq_change';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,300],'Name',figtitlestr{FigureIndex});
subplot(1,3,1)

r_seq_initial_sig =r_seq_initial(p_seq_initial <0.05);
r_seq_initial_nonsig =r_seq_initial(~(p_seq_initial <0.05));

edges = -1:0.1:1;

%Initial Learning
[N_r_sig, edges] = histcounts(r_seq_initial_sig, edges,'Normalization', 'probability');
[N_r_nonsig, edges] = histcounts(r_seq_initial_nonsig, edges,'Normalization', 'probability');

x_plot=edges(1:end-1)+0.05;
b = bar(x_plot,[N_r_sig',N_r_nonsig'],'stacked','LineWidth',3);
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
box off;

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Initial Learning');
ylabel('Proportion of sequences');
%xlabel('Correlation Coefficients');

subplot(1,3,2)

r_seq_20_sig =r_seq_20(p_seq_20 <0.05);
r_seq_20_nonsig =r_seq_20(~(p_seq_20 <0.05));

%edges = 0:0.1:1;

%After 20 days 
[N_r_sig, edges] = histcounts(r_seq_20_sig, edges,'Normalization', 'probability');
[N_r_nonsig, edges] = histcounts(r_seq_20_nonsig, edges,'Normalization', 'probability');

x_plot=edges(1:end-1)+0.05;
b = bar(x_plot,[N_r_sig',N_r_nonsig'],'stacked','LineWidth',3);
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
box off;

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('After 15 days');
ylabel('Proportion of sequences');
xlabel('Correlation Coefficients');


subplot(1,3,3)


plot(1:length(r_seq_mean),r_seq_mean,'-r','LineWidth',3);
hold on
errorbar(1:length(r_seq_mean),r_seq_mean,r_seq_sem,'-k','CapSize',0,'LineWidth',3)

segline = linspace(1,length(r_seq_mean),50);
plot(segline ,zeros(1,length(segline)),'--k','LineWidth',1.5);

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off
ylabel('Correlation Coefficients');
xlabel('Training days');


%% Comparison between fractal based sequence and scene based sequence 

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='CompareSceneFrac';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,300],'Name',figtitlestr{FigureIndex});


max_val = max(max(NTrialCri_Initial_frac),max(NTrialCri_Initial_scene));
edges = 0:10:max_val;

[N_frac, edges] = histcounts(NTrialCri_Initial_frac, edges,'Normalization', 'probability');
[N_scene, edges] = histcounts(NTrialCri_Initial_scene, edges,'Normalization', 'probability');


subplot(2,2,1)
x_plot=edges(1:end-1)+5;
b = bar(x_plot,N_frac,'LineWidth',3);
b.FaceColor=[0.5,0.5,0.5];

hold on
yrange = ylim;
plot(Mean_NCri_frac,yrange(2)+0.1,'v','markersize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(Mean_NCri_frac*ones(length(yrange(1):0.1:yrange(2)+0.1)),yrange(1):0.1:yrange(2)+0.1,'--k','LineWidth',3);



ylabel('Proportion of sequence');
title('Fractal based sequence');

box off;
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%xlabel('Number of trials to criterium');


subplot(2,2,3)
x_plot=edges(1:end-1)+5;
b = bar(x_plot,N_scene,'LineWidth',3);
b.FaceColor=[0.5,0.5,0.5];

hold on
yrange = ylim;
plot(Mean_NCri_scene,yrange(2)+0.1,'v','markersize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(Mean_NCri_scene*ones(length(yrange(1):0.1:yrange(2)+0.1)),yrange(1):0.1:yrange(2)+0.1,'--k','LineWidth',3);

ylabel('Proportion of sequence');
title('Scene based sequence');

box off

xlabel('Number of trials to criterium');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');




max_val = max(max(InitialLearning_CompleteTime_frac),max(InitialLearning_CompleteTime_scene ));
edges = 0:100:max_val;

[N_frac, edges] = histcounts(InitialLearning_CompleteTime_frac, edges,'Normalization', 'probability');
[N_scene, edges] = histcounts(InitialLearning_CompleteTime_scene,edges,'Normalization', 'probability');

subplot(2,2,2)
x_plot=edges(1:end-1)+50;
b = bar(x_plot,N_frac,'LineWidth',3);
b.FaceColor=[0.5,0.5,0.5];

hold on
yrange = ylim;
plot(Mean_CT_frac,yrange(2)+0.1,'v','markersize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(Mean_CT_frac*ones(length(yrange(1):0.1:yrange(2)+0.1)),yrange(1):0.1:yrange(2)+0.1,'--k','LineWidth',3);



ylabel('Proportion of sequence');
title('Fractal based sequence');

box off;
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%xlabel('Number of trials to criterium');


subplot(2,2,4)
x_plot=edges(1:end-1)+50;
b = bar(x_plot,N_scene,'LineWidth',3);
b.FaceColor=[0.5,0.5,0.5];

hold on
yrange = ylim;
plot(Mean_CT_scene,yrange(2)+0.1,'v','markersize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(Mean_CT_scene*ones(length(yrange(1):0.1:yrange(2)+0.1)),yrange(1):0.1:yrange(2)+0.1,'--k','LineWidth',3);

ylabel('Proportion of sequence');
title('Scene based sequence');

box off

xlabel('Hyperset completion time');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%% Comparison before gap and after gap
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='Comparison before gap and after gap';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,300],'Name',figtitlestr{FigureIndex});

subplot(2,2,1)
plot(NTrialCri_Initial_before_gap,NTrialCri_post_initial_gap,'ok','MarkerFaceColor','r','MarkerSize',10,'LineWidth',3);

max_val = max([NTrialCri_Initial_before_gap,NTrialCri_post_initial_gap]);
min_val = min([NTrialCri_Initial_before_gap,NTrialCri_post_initial_gap]);

hold on

plot([min_val:10:max_val],[min_val:10:max_val],'--k','LineWidth',3);
xlim([min_val-10,max_val+10]);
ylim([min_val-10,max_val+10]);

xlabel('Number of Trials to Criterium(Initial)');
ylabel('Number of Trials to Criterium(After Gap)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;

subplot(2,2,2)
plot(NTrialCri_last_before_gap,NTrialCri_post_initial_gap,'ok','MarkerFaceColor','r','MarkerSize',10,'LineWidth',3);
max_val = max([NTrialCri_last_before_gap,NTrialCri_post_initial_gap]);
min_val = min([NTrialCri_last_before_gap,NTrialCri_post_initial_gap]);
hold on
plot([min_val:10:max_val],[min_val:10:max_val],'--k','LineWidth',3);
xlim([min_val-10,max_val+10]);
ylim([min_val-10,max_val+10]);
box off

xlabel('Number of Trials to Criterium(Well Trained)');
ylabel('Number of Trials to Criterium(After Gap)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



subplot(2,2,3)
plot(c_time_initial_gap ,ctime_last_after_gap ,'ok','MarkerFaceColor','r','MarkerSize',10,'LineWidth',3);

max_val = max([c_time_initial_gap ,ctime_last_after_gap]);
min_val = min([c_time_initial_gap ,ctime_last_after_gap]);

hold on

plot([min_val:10:max_val],[min_val:10:max_val],'--k','LineWidth',3);
xlim([min_val-10,max_val+10]);
ylim([min_val-10,max_val+10]);

xlabel('Hyperset Completion Time(Initial)');
ylabel('Hyperset Completion Time(After Gap)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;

subplot(2,2,4)
plot(ctime_last_before_gap,ctime_last_after_gap,'ok','MarkerFaceColor','r','MarkerSize',10,'LineWidth',3);
max_val = max([ctime_last_before_gap,ctime_last_after_gap]);
min_val = min([ctime_last_before_gap,ctime_last_after_gap]);
hold on
plot([min_val:10:max_val],[min_val:10:max_val],'--k','LineWidth',3);
xlim([min_val-10,max_val+10]);
ylim([min_val-10,max_val+10]);
box off

xlabel('Hyperset Completion Time(Well Trained)');
ylabel('Hyperset Completion Time(After Gap)');


set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');







end %If plot figure
%Plot all the initial training curve 


%{
for i=1:length(FilesName)
    for j = 1:length(FilesName{i})
        clear OutputData;
        ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];
  
     load(FilesName{i}{j});
     
     OutputDataTemp=OutputData;
     Task(index)=OutputData.SequenceLearning.TaskCode;
     TaskName{index}=OutputData.SequenceLearning.Task;

     FilesNameAll{index}=FilesName_Each{i};
     
     
     Data =OutputData.SequenceLearning.DataStamp;

    

     Seq = Data('UniqueSeq');
     CorrectRate= Data('CorrectRate');

     TgtAcTime = Data('TgtAcTime');
     TotalCompleteTrial = Data('TotalCompleteTrial');

     CorrectTrial = Data('CorrectTrial');

     NTrialToCri{index} = Data('NTrialToCri');

     r_seq = Data('r_seq');

     p_seq{index} = Data('p_seq');

     TgtATimeEach_Mean{index} = Data('TgtATimeEach_Mean');

     

     for ss = 1:length(Seq)

        SeqString = string(Seq(ss));

        if CorrectTrial(ss)>3 && TotalCompleteTrial(ss)>5

            %Correct Rate

            if sum(contains(keys(seq_cr),SeqString))>0  

                seq_cr(string(Seq(ss)))=[seq_cr(string(Seq(ss))),CorrectRate(ss)];
            else
                seq_cr(string(Seq(ss)))=[CorrectRate(ss)];
            end

            %Target acquisition Time

            if sum(contains(keys(seq_tt),SeqString))>0  

                seq_tt(string(Seq(ss)))=[seq_tt(string(Seq(ss))),TgtAcTime(ss)];
            else
                seq_tt(string(Seq(ss)))=[TgtAcTime(ss)];
            end

            %Slope of TAT in each order
             if sum(contains(keys(seq_r),SeqString))>0  

                seq_r(string(Seq(ss)))=[seq_r(string(Seq(ss))),r_seq(ss)];
            else
                seq_r(string(Seq(ss)))=[r_seq(ss)];
             end


            


           



        end

     end








     index=index+1;




    end
     
     
end

SeqAll = str2double(keys(seq_cr));

PlotExample = 1;
Example = {'38','26','53'};


Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
Color=[Color;Color*0.5];
Color=[Color;Color*0.8];


if ShowFigureFlag
 
FigureStartNum=100;
FigureIndex=1;

% Learning curve
figtitlestr{FigureIndex}='PSTH_LearningCurve_CR';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

for i = 1:length(SeqAll)
    CR = seq_cr(string(SeqAll(i)));
    plot(1:length(CR),CR,'-ob','color',Color(i,:),'LineWidth',3);
    hold on



end
box off
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');
ylim([0,1]);
xlabel('Training Day');
ylabel('Correct Rate');
legend(string(SeqAll));

%Target acqusition time

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_LearningCurve_TargetAquistionTime';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

for i = 1:length(SeqAll)
    TT = seq_tt(string(SeqAll(i)));
    plot(1:length(TT),TT,'-ob','color',Color(i,:),'LineWidth',3);
    hold on



end
box off
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

xlabel('Training Day');
ylabel('Target Aqusition Time(ms)');
legend(string(SeqAll));

%TargetAqusitionTime Slope

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_LearningCurve_Slope';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

for i = 1:length(SeqAll)
    r = seq_r(string(SeqAll(i)));
    plot(1:length(r),r,'-ob','color',Color(i,:),'LineWidth',3);
    hold on



end
box off
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

xlabel('Training Day');
ylabel('Slope of the Target Aquisition Time(ms)');
legend(string(SeqAll));



if PlotExample
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_CR';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

    for i = 1:length(Example)
    plot(1:length(seq_cr(Example{i})),seq_cr(Example{i}),'-ob','color',Color(i,:),'LineWidth',3);
    hold on

    box off
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');
    ylim([0,1]);
    xlabel('Training Day');
    ylabel('Correct Rate');
    legend(Example);


    end

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_TT';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

    for i = 1:length(Example)
    plot(1:length(seq_tt(Example{i})),seq_tt(Example{i}),'-ob','color',Color(i,:),'LineWidth',3);
    hold on

    box off
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');
    
    xlabel('Training Day');
    ylabel('Target Aqusition Time(ms)');
    legend(Example);


    end

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_r';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

    for i = 1:length(Example)
    plot(1:length(seq_r(Example{i})),seq_r(Example{i}),'-ob','color',Color(i,:),'LineWidth',3);
    hold on

    box off
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');
    
    xlabel('Training Day');
    ylabel('Slope of the target acqusition time(ms)');
    legend(Example);


    end


end
%}

%{
figure
for i = 1:length(SeqAll)
    CR = seq_cr(string(SeqAll(i)));
    plot(CR,'LineWidth',3);
    hold on



end
legend(string(SeqAll));
%}

end
function seq = fillInNaN(c);

maxnum = max(cellfun(@numel,c));

seq = cell2mat(cellfun(@(x) [x;NaN*ones(maxnum-length(x),1)],c,'UniformOutput',0))';



end