function [X, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
    %KALMANFILTER Filters measurements sequence Y using a Kalman filter.
    %
    %Input:
    %   Y           [m x N] Measurement sequence
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   A           [n x n] State transition matrix
    %   Q           [n x n] Process noise covariance
    %   H           [m x n] Measurement model matrix
    %   R           [m x m] Measurement noise covariance
    %
    %Output:
    %   x           [n x N] Estimated state vector sequence
    %   P           [n x n x N] Filter error convariance
    %

    %% Parameters
    N = size(Y,2);

    n = length(x_0);
    m = size(Y,1);

    %% Data allocation
    X = zeros(n,N);
    P = zeros(n,n,N);

    [x_pred, p_pred] = linearPrediction(x_0, P_0, A, Q);
    [X(:,1), P(:,:,1)] = linearUpdate(x_pred, p_pred, Y(:,1), H, R);

    for i = 2:N
        [x_pred, p_pred] = linearPrediction(X(:,i-1), P(:,:,i-1), A, Q);
        [X(:,i), P(:,:,i)] = linearUpdate(x_pred, p_pred, Y(:,i), H, R);
    end

    end

