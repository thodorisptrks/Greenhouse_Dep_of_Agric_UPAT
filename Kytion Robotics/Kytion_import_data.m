%% data concatenation for KYTION
clear all

cd 'G:\Το Drive μου\GitHub\Greenhouse_UPAT\Kytion Robotics'

projectdir = 'G:\Το Drive μου\Data_Greenhouse\KYTION';
dinfo = dir(fullfile(projectdir, '*.txt'));   %use appropriate extension
kytion_data = readtable((fullfile(projectdir,dinfo.name)),'Delimiter',{',','"'},'ReadVariableNames', false);

kytion_data = kytion_data(:,["Var4","Var7","Var12","Var15"]);

kytion_data.TimeStamp = kytion_data.Var4 + kytion_data.Var7;
kytion_data.TimeStamp.Format = 'dd/MM/yyyy HH:mm:ss';

kytion_data = table(kytion_data.TimeStamp,kytion_data.Var12,kytion_data.Var15);
kytion_data.Properties.VariableNames = {'TimeStamp','Temperature','Rel_Humidity'};

starting_datetime_idx = find(kytion_data.TimeStamp == '12/09/2022 12:14:49');
kytion_data = kytion_data(starting_datetime_idx:end,:);

kytion_data = table2timetable(kytion_data);

kytion_data = retime(kytion_data,'regular','mean','TimeStep',minutes(10));
