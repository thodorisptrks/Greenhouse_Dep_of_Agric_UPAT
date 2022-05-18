clear all
clc
close all

cd 'G:\Το Drive μου\Scripts'

%% import data
data = readtable('data_stations.csv');

%% subsets indices

% training subset indices before testing subset
training_data_start_idx = find(data.TIMESTAMP == data.TIMESTAMP(1));
training_data_pause_idx = find(data.TIMESTAMP == datetime('2022-03-27 6:10'));   

training_idx_range_before = [training_data_start_idx:training_data_pause_idx];

% testing subset indices
test_data_start_idx = find(data.TIMESTAMP == datetime('2022-03-27 6:20'));
test_data_end_idx = find(data.TIMESTAMP == datetime('2022-03-30 6:20'));

testing_idx_range = [test_data_start_idx:test_data_end_idx];

% training subset indices after testing subset
training_data_cont_idx = find(data.TIMESTAMP == datetime('2022-03-30 6:30'));
training_data_end_idx = find(data.TIMESTAMP == data.TIMESTAMP(end));

training_idx_range_after = [training_data_cont_idx:training_data_end_idx];

%% Figures

figure1 = figure('WindowState','maximized','Color',[1 1 1]);
plot(training_idx_range_before,data.Tair_Avg(training_idx_range_before),'-',...
    'LineWidth',1.7,'Color','b');
hold on
plot(testing_idx_range, data.Tair_Avg(testing_idx_range),'-',...
    'LineWidth',1.7,'Color','r');
hold on
plot(training_idx_range_after,data.Tair_Avg(training_idx_range_after),'-',...
    'LineWidth',1.7,'Color','b');
hold off
xlabel('Number of Samples','FontWeight','bold')
xlim([1,training_data_end_idx]);
xticks(1:500:training_data_end_idx)
ylabel('Outside Temperature (^oC)','FontWeight','bold')
ylim([0 30])
yticks(0:5:30)
legend('Training/Validation data','Test data')
ax = gca;
ax.FontSize = 16; 
grid on
box on
saveas(gca,'G:\Το Drive μου\Scripts\graphs\T_out.png') % save timetable as csv

figure2 = figure('WindowState','maximized','Color',[1 1 1]);
plot(training_idx_range_before,data.RH_Avg(training_idx_range_before),'-',...
    'LineWidth',1.7,'Color','b');
hold on
plot(testing_idx_range, data.RH_Avg(testing_idx_range),'-',...
    'LineWidth',1.7,'Color','r');
hold on
plot(training_idx_range_after,data.RH_Avg(training_idx_range_after),'-',...
    'LineWidth',1.7,'Color','b');
hold off
xlabel('Number of Samples','FontWeight','bold')
xlim([1,training_data_end_idx]);
xticks(1:500:training_data_end_idx)
ylabel('Outside Relative Humidity (%)','FontWeight','bold')
ylim([10 110])
yticks(10:10:110)
legend('Training/Validation data','Test data')
ax = gca;
ax.FontSize = 16; 
grid on
box on
saveas(gca,'G:\Το Drive μου\Scripts\graphs\RH_out.png') % save timetable as csv

figure3 = figure('WindowState','maximized','Color',[1 1 1]);
plot(training_idx_range_before,data.WS_Avg(training_idx_range_before),'-',...
    'LineWidth',1.7,'Color','b');
hold on
plot(testing_idx_range, data.WS_Avg(testing_idx_range),'-',...
    'LineWidth',1.7,'Color','r');
hold on
plot(training_idx_range_after,data.WS_Avg(training_idx_range_after),'-',...
    'LineWidth',1.7,'Color','b');
hold off
xlabel('Number of Samples','FontWeight','bold')
xlim([1,training_data_end_idx]);
xticks(1:500:training_data_end_idx)
ylabel('Wind Speed (m*s^-^1)','FontWeight','bold')
ylim([0 7.5])
yticks(0:1:7)
legend('Training/Validation data','Test data')
ax = gca;
ax.FontSize = 16; 
grid on
box on
saveas(gca,'G:\Το Drive μου\Scripts\graphs\WS.png') % save as png

figure4 = figure('WindowState','maximized','Color',[1 1 1]);
plot(training_idx_range_before,data.SolRad_Out_Avg(training_idx_range_before),'-',...
    'LineWidth',1.7,'Color','b');
hold on
plot(testing_idx_range, data.SolRad_Out_Avg(testing_idx_range),'-',...
    'LineWidth',1.7,'Color','r');
hold on
plot(training_idx_range_after,data.SolRad_Out_Avg(training_idx_range_after),'-',...
    'LineWidth',1.7,'Color','b');
hold off
xlabel('Number of Samples','FontWeight','bold')
xlim([1,training_data_end_idx]);
xticks(1:500:training_data_end_idx)
ylabel('Solar Irradiance (W*m^-^2)','FontWeight','bold')
ylim([0 1200])
yticks(0:100:1200)
legend('Training/Validation data','Test data')
ax = gca;
ax.FontSize = 16; 
grid on
box on
saveas(gca,'G:\Το Drive μου\Scripts\graphs\Solar_Irrad.png') % save as png

figure5 = figure('WindowState','maximized','Color',[1 1 1]);
plot(training_idx_range_before,data.Tin_Avg(training_idx_range_before),'-',...
    'LineWidth',1.7,'Color','b');
hold on
plot(testing_idx_range, data.Tin_Avg(testing_idx_range),'-',...
    'LineWidth',1.7,'Color','r');
hold on
plot(training_idx_range_after,data.Tin_Avg(training_idx_range_after),'-',...
    'LineWidth',1.7,'Color','b');
hold off
xlabel('Number of Samples','FontWeight','bold')
xlim([1,training_data_end_idx]);
xticks(1:500:training_data_end_idx)
ylabel('Indoor Temperature (^oC)','FontWeight','bold')
ylim([0 60])
yticks(0:5:60)
legend('Training/Validation data','Test data')
ax = gca;
ax.FontSize = 16; 
grid on
box on
saveas(gca,'G:\Το Drive μου\Scripts\graphs\T_in.png') % save as png

figure6 = figure('WindowState','maximized','Color',[1 1 1]);
plot(training_idx_range_before,data.RHin_Avg(training_idx_range_before),'-',...
    'LineWidth',1.7,'Color','b');
hold on
plot(testing_idx_range, data.RHin_Avg(testing_idx_range),'-',...
    'LineWidth',1.7,'Color','r');
hold on
plot(training_idx_range_after,data.RHin_Avg(training_idx_range_after),'-',...
    'LineWidth',1.7,'Color','b');
hold off
xlabel('Number of Samples','FontWeight','bold')
xlim([1,training_data_end_idx]);
xticks(1:500:training_data_end_idx)
ylabel('Indoor Relative Humidity (%)','FontWeight','bold')
ylim([0 100])
yticks(0:10:100)
legend('Training/Validation data','Test data')
ax = gca;
ax.FontSize = 16; 
grid on
box on
saveas(gca,'G:\Το Drive μου\Scripts\graphs\RH_in.png') % save as png

