function [max_error_T,max_error_RH,mae_T,mae_RH,rmse_T,rmse_RH,r2_T,r2_RH,errors_T,errors_RH] =...
    NN_thodoris(hid_nodes,trn_alg,trnsf_fcn_1,trnsf_fcn_2,tr_ratio,val_ratio,tr_input,tr_target,tst_input,Temp_obs,RH_obs,...
    rmax,rmin,temp_max,temp_min,rh_max,rh_min)

    rng('default')
    hiddenLayerSize = hid_nodes;
    trainFcn = trn_alg;   % training algorithm 
    net = feedforwardnet(hiddenLayerSize,trainFcn);   % Feedforward NN
    net.layers{1}.transferFcn = trnsf_fcn_1;   % activation fcn in hidden
    net.layers{2}.transferFcn = trnsf_fcn_2;   % activation fcn in output
    net.divideFcn = 'dividerand';   % random separation
    net.divideParam.trainRatio = tr_ratio;   % training ratio
    net.divideParam.valRatio = val_ratio;   % validation ratio
    net.divideParam.testRatio = 0;
    net.performFcn = 'mse';   % perform fcn

    net = configure(net,tr_input,tr_target);   % configure the network

    [net,tr] = train(net,tr_input,tr_target);   % train the network

    % test the network
    y_test_predictions = net(tst_input);   % predictions

    % Temperature
    Temp_t_pred = ((((y_test_predictions(1,:))' - rmin)./(rmax - rmin)).*(temp_max...
        - temp_min)) + temp_min;   % de-normalize the temperature predictions
    
    Temp_t_obs = Temp_obs;   % observations
    
    % statistics for each struct
    
    abs_errors_T = abs(Temp_t_pred - Temp_t_obs);   % absolute errors between predictions & observations
    max_error_T = max(abs_errors_T,[],1);   % max absolute error
    mae_T = mean(abs(Temp_t_pred - Temp_t_obs)); % mean absolute error between predictions & observations
    rmse_T = sqrt(mean((Temp_t_pred - Temp_t_obs).^2));   % root mean square error between predictions & observations
    r2_T = 1 - ((sum((Temp_t_obs - Temp_t_pred).^2))/(sum((Temp_t_obs - mean(Temp_t_obs)).^2))); % r-squared
    errors_T = Temp_t_pred - Temp_t_obs;    % errors between predictions & observations
    
    % Relative Humidity
    RH_t_pred = ((((y_test_predictions(2,:))' - rmin)./(rmax - rmin)).*(rh_max...
        - rh_min)) + rh_min;   % de-normalize the rel. humidity predictions
    
    RH_t_obs = RH_obs;   % observations
    
    abs_errors_RH = abs(RH_t_pred - RH_t_obs);   % absolute errors between predictions & observations
    max_error_RH = max(abs_errors_RH,[],1);   % max absolute error
    mae_RH = mean(abs(RH_t_pred - RH_t_obs));   % mean absolute error between predictions & observations
    rmse_RH = sqrt(mean((RH_t_pred - RH_t_obs).^2));    % root mean square error between predictions & observations
    r2_RH = 1 - ((sum((RH_t_obs - RH_t_pred).^2))/(sum((RH_t_obs - mean(RH_t_obs)).^2)));   % r-squared
    errors_RH = Temp_t_pred - Temp_t_obs;   % errors between predictions & observations
end