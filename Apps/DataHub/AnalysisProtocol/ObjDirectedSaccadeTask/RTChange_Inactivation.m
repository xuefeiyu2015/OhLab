function RTChange_Inactivation(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag)
%Export to local file
%keyboard
%{
export_path = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Data/FEF_Inactivation_OHLab';
% Prepare data for export
InjectionCondition = zeros(Data.TrialNum,1);%Before injection
InjectionCondition(381:540) = NaN; %During injection
InjectionCondition(541:1341) = 1;%After injection
InjectionCondition(1341:end) = 2;%After 2 injection

Data.InjectionCondition = InjectionCondition;

output_file = fullfile(export_path, 'Adams102125FEFInactivation_FRAC.mat');
save(output_file, 'Data');
%disp('file exported')
%}

end