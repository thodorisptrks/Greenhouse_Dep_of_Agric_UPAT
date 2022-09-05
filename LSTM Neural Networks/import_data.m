clear all
clc
close all

%% import data from .txt files

projectdir = 'G:\Το Drive μου\Data_Greenhouse\txt_files';
dinfo = dir(fullfile(projectdir, '*.txt'));   %use appropriate extension
filenames = fullfile(projectdir, {dinfo.name});
nfiles = length(filenames);
tables = cell(nfiles,1);

for K = 1 : nfiles
    tables{K} = readtable(filenames{K}); % create tables from files
end

data_stations = vertcat(tables{:}); % put them all in a single table
data_stations = rmmissing(data_stations); %remove rows with NaN
data_stations = table2timetable(data_stations); % convert to timetable

dt = minutes(10); % set the duration
data_stations = retime(data_stations,'regular','fillwithmissing','TimeStep',dt); % fill missing dates with NaN

data_stations = data_stations(:,["Tair_Avg","RH_Avg","WS_Avg","SolRad_Out_Avg",...
    "Tair_in_1_Avg","RH_in_1_Avg","Sol_in_1_Avg","Par_in_1_Avg","Tair_in_2_Avg","RH_in_2_Avg",...
    "Sol_in_2_Avg","Par_in_2_Avg"]); % keep needed variables

%% average of inside temperature and relative humidity

data_stations.Tin_Avg = mean([data_stations.Tair_in_1_Avg,data_stations.Tair_in_2_Avg],2);
data_stations.RHin_Avg = mean([data_stations.RH_in_1_Avg,data_stations.RH_in_2_Avg],2);
data_stations.Solin_Avg = mean([data_stations.Sol_in_1_Avg,data_stations.Sol_in_2_Avg],2);
data_stations.Par_Avg = mean([data_stations.Par_in_1_Avg,data_stations.Par_in_2_Avg],2);

%% save dataset as csv 

writetimetable(data_stations,'G:\Το Drive μου\GitHub\Greenhouse_UPAT\LSTM Neural Networks\data_stations.csv','Delimiter',';')
