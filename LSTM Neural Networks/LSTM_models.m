clear all
clc
close all

cd 'G:\Το Drive μου\GitHub\Greenhouse_UPAT\LSTM Neural Networks'

%%
dataset = readtimetable("data_stations.csv");
dataset = rmmissing(dataset);

data_train = dataset(1:end-1008,:);
data_test = dataset(end-1007:end,:);

data_train = [data_train.Tair_Avg,data_train.RH_Avg,data_train.WS_Avg,data_train.SolRad_Out_Avg,((data_train.Tair_in_1_Avg + data_train.Tair_in_2_Avg)./2)];
data_test = [data_test.Tair_Avg,data_test.RH_Avg,data_test.WS_Avg,data_test.SolRad_Out_Avg,((data_test.Tair_in_1_Avg + data_test.Tair_in_2_Avg)./2)];

%% normalization according to mean and std after the partition to train and test

mu_train = mean(data_train,1);
std_train = std(data_train,1);

mu_test = mean(data_test,1);
std_test = std(data_test,1);

data_train_norm = (data_train - mu_train)./std_train;
data_test_norm = (data_test - mu_test)./std_test;

inputs_train = num2cell((data_train_norm(:,1:4))',1);
inputs_test = num2cell((data_test_norm(:,1:4))',1);

outputs_train = num2cell((data_train_norm(:,5))',1);
outputs_test = num2cell((data_test_norm(:,5))',1);

%% Define Network Architecture
featureDimension = size(inputs_train{1},1);
numResponses = size(outputs_train{1},1);
numHiddenUnits = 500;

layers = [ ...
    sequenceInputLayer(featureDimension)
    lstmLayer(numHiddenUnits,'OutputMode','sequence')
    fullyConnectedLayer(numResponses)
    regressionLayer];

maxepochs = 100;
% miniBatchSize = 100;

options = trainingOptions('adam', ...  %%adam
    'MaxEpochs',maxepochs, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',0, ...
    'Plots','training-progress');

%% Train the Network
net = trainNetwork(inputs_train,outputs_train,layers,options);

%% Test the Network

net = resetState(net);
Predicts = predict(net,inputs_test);

Predicts = (cell2mat(Predicts).*std_test(1,5)) + mu_test(1,5);
Observs = (cell2mat(outputs_test').*std_test(1,5)) + mu_test(1,5);

abs_errors = abs(Predicts - Observs);   % absolute errors between predictions & observations
max_error = max(abs_errors,[],1);   % max absolute error
mae = mean(abs(Predicts - Observs));   % mean absolute error between predictions & observations
rmse = sqrt(mean((Predicts - Observs).^2));    % root mean square error between predictions & observations
r2 = 1 - ((sum((Observs - Predicts).^2))/(sum((Observs - mean(Observs)).^2)));   % r-squared
errors = Predicts - Observs;   % errors between predictions & observations

% 
figure1 = figure('WindowState','maximized','Color','w');
plot(Observs,'-*','Color','b')
hold on
plot(Predicts,'-o','Color','r')
ylabel('Indoor Temperature')
xlabel('Number of Samples')
xlim([0 1007])
legend('Observations','Predictions','Location','best')
grid on 
box on

figure2 = figure('WindowState','maximized','Color','w');
plot(Observs,Predicts,'ob','MarkerFaceColor','b')
ylabel('Predictions')
xlabel('Observations')
r = refline(1,0);
r.LineWidth = 2;
r.Color = 'r';
r.LineStyle = '--';
grid on 
box on

figure3 = figure('WindowState','maximized','Color','w');
stem(errors)
ylabel('Errors (Pred. - Obs.)','FontWeight','bold')
xlim([0 1007])
xlabel('Number of Samples','FontWeight','bold')
grid on 
box on

figure4 = figure('WindowState','maximized','Color','w');
yyaxis left
stem(errors)
ylabel('Errors (Pred. - Obs.)','FontWeight','bold')
xlim([0 1007])
xlabel('Number of Samples','FontWeight','bold')
grid on 
box on
yyaxis right
plot(Observs,'-','Color','b')
hold on
plot(Predicts,'-','Color','r')
ylabel('Indoor Temperature')
legend('Observations','Predictions')
