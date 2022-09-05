clear all
clc
close all

cd 'G:\Το Drive μου\Scripts'

%% import data
data = readtable('data_stations.csv');

%% delay for indoor temperature and relative humidity

data.T_in_Avg_10min_back = [data.Tin_Avg(2:end,1);NaN];   % one timestep before
data.T_in_Avg_20min_back = [data.Tin_Avg(3:end,1);NaN;NaN];   % two timesteps before 
data.T_in_Avg_30min_back = [data.Tin_Avg(4:end,1);NaN;NaN;NaN];   % three timesteps before 

data.RH_in_Avg_10min_back = [data.RHin_Avg(2:end,1);NaN];   % one timestep before
data.RH_in_Avg_20min_back = [data.RHin_Avg(3:end,1);NaN;NaN];   % two timesteps before 
data.RH_in_Avg_30min_back = [data.RHin_Avg(4:end,1);NaN;NaN;NaN];   % three timesteps before 

%%
sol_less_4_5 = data.SolRad_Out_Avg < 4.5;
data.SolRad_Out_Avg(sol_less_4_5) = 0;

%%
data = rmmissing(data); % remove rows with NaN

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

%%
data = data(:,2:end);   % delete TIMESTAMP

data_training = data([training_idx_range_before,training_idx_range_after],:);   % training dataset
data_testing = data(testing_idx_range,:);   % testing dataset

data_training_array = table2array(data_training);   % turn table to array
data_testing_array = table2array(data_testing);   % turn table to array

%%  normalization of training and testing dataset in a range of 0.01 to 0.99

[data_training_mins,data_training_maxes,data_training_norm] = normlzn(data_training_array,0.99,0.01);

[data_testing_mins,data_testing_maxes,data_testing_norm] = normlzn(data_testing_array,0.99,0.01);

%% Inputs

% training
I_training = data_training_array(:,[1 2 3 4 10 11 12 13 14 15]);
I_training_norm = data_training_norm(:,[1 2 3 4 10 11 12 13 14 15]);   % normalized
I_training_norm_tr = I_training_norm';   % transposed

% testing
I_testing = data_testing_array(:,[1 2 3 4 10 11 12 13 14 15]);
I_testing_norm = data_testing_norm(:,[1 2 3 4 10 11 12 13 14 15]);   % normalized
I_testing_norm_tr = I_testing_norm';   % transposed

%% Targets 

% training
T_training = data_training_array(:,[6 7]);
T_training_norm = data_training_norm(:,[6 7]);   % normalized
T_training_norm_tr = T_training_norm';   % transposed

% testing
T_testing = data_testing_array(:,[6 7]);
T_testing_norm = data_testing_norm(:,[6 7]);   % normalized
T_testing_norm_tr = T_testing_norm';   % transposed

%% Inputs for the function .... mlp_nn(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
    % 1. number of hidden nodes
    % 2. training algorithm
    % 3. transfer function for hidden layer
    % 4. transfer function for output layer
    % 5. training ratio
    % 6. validation ratio
    % 7. inputs for training
    % 8. targets for training
    % 9. inputs for testing
    % 10. Temperature observations
    % 11. Relative Humidity observations
    % 12. rmax for de-normalization
    % 13. rmin for de-normalization
    % 14. Temperature max for de-normalization
    % 15. Temperature min for de-normalization
    % 16. Relative Humidity max for de-normalization
    % 17. Relative Humidity min for de-normalization

%% Outputs from the function .... (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
    % 1. Temperature MAX ERROR
    % 2. Relative Humidity MAX ERROR
    % 3. Temperature MEAN ABSOLUTE ERROR
    % 4. Relative Humidity MEAN ABSOLUTE ERROR
    % 5. Temperature ROOT MEAN SQUARE ERROR
    % 6. Relative Humidity ROOT MEAN SQUARE ERROR
    % 7. Temperature R-SQUARED
    % 8. Relative Humidity R-SQUARED
    % 9. Temperature ERRORS
    % 10. Relative Humidity ERRORS

%% activation fcn in hidden --> logistic sigmoid
% training algorithm --> Levenberg-Marquardt backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainlm','logsig','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> Bayesian Regularization backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbr','logsig','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> BFGS Quasi Newton backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbfg','logsig','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

%% activation fcn in hidden --> radial basis
% training algorithm --> Levenberg-Marquardt backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainlm','radbas','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> Bayesian Regularization backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbr','radbas','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> BFGS Quasi Newton backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbfg','radbas','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

%% activation fcn in hidden --> positive linear
% training algorithm --> Levenberg-Marquardt backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainlm','poslin','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> Bayesian Regularization backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbr','poslin','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> BFGS Quasi Newton backpropagation
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbfg','poslin','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,0.01,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

%%  normalization of training and testing dataset in a range of -0.99 to 0.99

[data_training_mins,data_training_maxes,data_training_norm] = normlzn(data_training_array,0.99,-0.99);

[data_testing_mins,data_testing_maxes,data_testing_norm] = normlzn(data_testing_array,0.99,-0.99);

%% Inputs

% training
I_training = data_training_array(:,[1 2 3 4 10 11 12 13 14 15]);
I_training_norm = data_training_norm(:,[1 2 3 4 10 11 12 13 14 15]);   % normalized
I_training_norm_tr = I_training_norm';   % transposed

% testing
I_testing = data_testing_array(:,[1 2 3 4 10 11 12 13 14 15]);
I_testing_norm = data_testing_norm(:,[1 2 3 4 10 11 12 13 14 15]);   % normalized
I_testing_norm_tr = I_testing_norm';   % transposed

%% Targets 

% training
T_training = data_training_array(:,[6 7]);
T_training_norm = data_training_norm(:,[6 7]);   % normalized
T_training_norm_tr = T_training_norm';   % transposed

% testing
T_testing = data_testing_array(:,[6 7]);
T_testing_norm = data_testing_norm(:,[6 7]);   % normalized
T_testing_norm_tr = T_testing_norm';   % transposed

%%
% training algorithm --> Levenberg-Marquardt backpropagation
% activation fcn in hidden --> hyperbolic tangent
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainlm','tansig','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,-0.99,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> Bayesian Regularization backpropagation
% activation fcn in hidden --> hyperbolic tangent
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbr','tansig','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,-0.99,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end

% training algorithm --> BFGS Quasi Newton backpropagation
% activation fcn in hidden --> hyperbolic tangent
% activation fcn in output --> linear
% nodes in hidden --> 1 to 20 
% training ratio --> 80%
% validation ratio --> 20%
% testing for the testing subset of dataset
% perform fcn --> mean square error

for i=1:20
    [max_error_T(:,i),max_error_RH(:,i),mae_T(:,i),mae_RH(:,i),rmse_T(:,i),rmse_RH(:,i),...
        r2_T(:,i),r2_RH(:,i),errors_T(:,i),errors_RH(:,i)] = ...
        mlp_nn(i,'trainbfg','tansig','purelin',0.8,0.2,I_training_norm_tr,T_training_norm_tr,...
        I_testing_norm_tr,T_testing(:,1),T_testing(:,2),0.99,-0.99,data_testing_maxes(1,6),...
        data_testing_mins(1,6),data_testing_maxes(1,7),data_testing_mins(1,7));
end
