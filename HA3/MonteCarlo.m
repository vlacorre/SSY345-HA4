function [x_err, P_err] = MonteCarlo(x_0,P_0,f,Q,h,R,type,MC,N)
    % Allocate memory
    x_err = zeros(length(x_0),MC);
    P_err = zeros(length(x_0),MC);
    
    for imc = 1:MC
        % Simulate state sequence
        X = genNonLinearStateSequence(x_0, P_0, f, Q, N);
        % Simulate measurements
        Y = genNonLinearMeasurementSequence(X, h, R);
        % Run Kalman filter (you need to run all three, for comparison)
        [xf,Pf,xp,Pp] = nonLinearKalmanFilter(Y,x_0,P_0,f,Q,h,R,type);
        % Save the estimation errors and the prediction errors!
        x_err(:,imc) = sqrt(sum((X(:,2:end)-xf).^2,2));
        P_err(:,imc) = sqrt(sum((X(:,2:end)-xp).^2,2));
    end
end