function [xs, Ps, xf, Pf, xp, Pp] = ...
    nonLinRTSsmoother(Y, x_0, P_0, f, T, Q, S, h, R, sigmaPoints, type)
    %NONLINRTSSMOOTHER Filters measurement sequence Y using a
    % non-linear Kalman filter.
    %
    %Input:
    %   Y           [m x N] Measurement sequence for times 1,...,N
    %   x_0         [n x 1] Prior mean for time 0
    %   P_0         [n x n] Prior covariance
    %   f                   Motion model function handle
    %   T                   Sampling time
    %   Q           [n x n] Process noise covariance
    %   S           [n x N] Sensor position vector sequence
    %   h                   Measurement model function handle
    %   R           [n x n] Measurement noise covariance
    %   sigmaPoints Handle to function that generates sigma points.
    %   type        String that specifies type of non-linear filter/smoother
    %
    %Output:
    %   xf          [n x N]     Filtered estimates for times 1,...,N
    %   Pf          [n x n x N] Filter error convariance
    %   xp          [n x N]     Predicted estimates for times 1,...,N
    %   Pp          [n x n x N] Filter error convariance
    %   xs          [n x N]     Smoothed estimates for times 1,...,N
    %   Ps          [n x n x N] Smoothing error convariance

    %% Parameters
    N = size(Y, 2);
    n = size(x_0, 1);

    % Filter
    xf = zeros(n,N);
    Pf = zeros(n,n,N);
    xp = zeros(n,N);
    Pp = zeros(n,n,N);
    [xp(:,1), Pp(:,:,1)] = nonLinKFprediction(x_0, P_0, f, T, Q, sigmaPoints, type);
    [xf(:,1), Pf(:,:,1)] = nonLinKFupdate(xp(:,1), Pp(:,:,1), Y(:,1), S, h, R, sigmaPoints, type);

    for i = 2:N
        [xp(:,i), Pp(:,:,i)] = nonLinKFprediction(xf(:,i-1), Pf(:,:,i-1), f, T, Q, sigmaPoints, type);
        [xf(:,i), Pf(:,:,i)] = nonLinKFupdate(xp(:,i), Pp(:,:,i), Y(:,i), S, h, R, sigmaPoints, type);
    end


    % Smoother
    xs = zeros(n,N);
    Ps = zeros(n,n,N);

    xs(:,N) = xf(:,N);
    Ps(:,:,N) = Pf(:,:,N);

    for i = N-1:-1:1
        [xs(:,i), Ps(:,:,i)] = nonLinRTSSupdate(xs(:,i+1), Ps(:,:,i+1), xf(:,i), Pf(:,:,i), xp(:,i+1), Pp(:,:,i+1), f, T, sigmaPoints, type);
    end

end





% GIVEN FUNCTIONS

function [xs, Ps] = nonLinRTSSupdate(xs_kplus1, ...
    Ps_kplus1, ...
    xf_k, ...
    Pf_k, ...
    xp_kplus1, ...
    Pp_kplus1, ...
    f, ...
    T, ...
    sigmaPoints, ...
    type)
    %NONLINRTSSUPDATE Calculates mean and covariance of smoothed state
    % density, using a non-linear Gaussian model.
    %
    %Input:
    %   xs_kplus1   Smooting estimate for state at time k+1
    %   Ps_kplus1   Smoothing error covariance for state at time k+1
    %   xf_k        Filter estimate for state at time k
    %   Pf_k        Filter error covariance for state at time k
    %   xp_kplus1   Prediction estimate for state at time k+1
    %   Pp_kplus1   Prediction error covariance for state at time k+1
    %   f           Motion model function handle
    %   T           Sampling time
    %   sigmaPoints Handle to function that generates sigma points.
    %   type        String that specifies type of non-linear filter/smoother
    %
    %Output:
    %   xs          Smoothed estimate of state at time k
    %   Ps          Smoothed error convariance for state at time k

    n = size(xs_kplus1,1);

    switch type
        case 'EKF'
            [x, Jac] = f(xf_k, T);
            Pk_kplus1_k = Pf_k*Jac';

        case 'UKF'
            [Chis, Ws] = sigmaPoints(xf_k, Pf_k, 'UKF');

            Pk_kplus1_k = 0;
            for i=1:2*n+1
                Pk_kplus1_k = Pk_kplus1_k + (Chis(:,i) - xf_k)*(f(Chis(:,i), T) - xp_kplus1)'*Ws(i);
            end

            % Make sure the cross-covariance matrix is semi-definite
            if min(eig(Pk_kplus1_k))<=0
                [v,e] = eig(Pk_kplus1_k, 'vector');
                e(e<0) = 1e-4;
                Pk_kplus1_k = v*diag(e)/v;
            end

        case 'CKF'
            [Chis, Ws] = sigmaPoints(xf_k, Pf_k, 'CKF');

            Pk_kplus1_k = 0;
            for i=1:2*n
                Pk_kplus1_k = Pk_kplus1_k + (Chis(:,i) - xf_k)*(f(Chis(:,i), T) - xp_kplus1)'*Ws(i);
            end

        otherwise
        error('Incorrect type of non-linear RTS smoothing')
    end

    % RTS aglorithm
    % Gk = Pk_k*Ak'/Pp_kplus1; % Linear formula for Gk
    Gk = Pk_kplus1_k/Pp_kplus1; % Non linear formula for Gk

    xs = xf_k + Gk*(xs_kplus1 - xp_kplus1);
    Ps = Pf_k - Gk*(Pp_kplus1 - Ps_kplus1)*Gk';

end

function [x, P] = nonLinKFprediction(x, P, f, T, Q, sigmaPoints, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%   T           Sampling time
%   Q           [n x n] Process noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

switch type
case 'EKF'

% Evaluate motion model
[fx, Fx] = f(x,T);
% State prediction
x = fx;
% Covariance prediciton
P = Fx*P*Fx' + Q;
% Make sure P is symmetric
P = 0.5*(P + P');

case 'UKF'

% Predict
[x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type);

if min(eig(P))<=0
[v,e] = eig(P);
emin = 1e-3;
e = diag(max(diag(e),emin));
P = v*e*v';
end

case 'CKF'

% Predict
[x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type);

otherwise
error('Incorrect type of non-linear Kalman filter')
end
end

function [x, P] = nonLinKFupdate(x, P, y, s, h, R, sigmaPoints, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   s           [2 x 1] sensor position vector
%   h           Measurement model function handle
%   R           [n x n] Measurement noise covariance
%   sigmaPoints Handle to function that generates sigma points.
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%


switch type
case 'EKF'

% Evaluate measurement model
[hx, Hx] = h(x,s);

% Innovation covariance
S = Hx*P*Hx' + R;
% Kalman gain
K = (P*Hx')/S;

% State update
x = x + K*(y - hx);
% Covariance update
P = P - K*S*K';

% Make sure P is symmetric
P = 0.5*(P + P');

case 'UKF'

% Update mean and covariance
[x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s, h, R, sigmaPoints, type);

if min(eig(P))<=0
[v,e] = eig(P);
emin = 1e-3;
e = diag(max(diag(e),emin));
P = v*e*v';
end

case 'CKF'

% Update mean and covariance
[x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s, h, R, sigmaPoints, type);

otherwise
error('Incorrect type of non-linear Kalman filter')
end

end


function [x, P] = predictMeanAndCovWithSigmaPoints(x, P, f, T, Q, sigmaPoints, type)
%
%PREDICTMEANANDCOVWITHSIGMAPOINTS computes the predicted mean and covariance
%
%Input:
%   x           [n x 1] mean vector
%   P           [n x n] covariance matrix
%   f           measurement model function handle
%   T           sample time
%   Q           [m x m] process noise covariance matrix
%
%Output:
%   x           [n x 1] Updated mean
%   P           [n x n] Updated covariance
%

% Compute sigma points
[SP,W] = sigmaPoints(x, P, type);

% Dimension of state and number of sigma points
[n, N] = size(SP);

% Allocate memory
fSP = zeros(n,N);

% Predict sigma points
for i = 1:N
[fSP(:,i),~] = f(SP(:,i),T);
end

% Compute the predicted mean
x = sum(fSP.*repmat(W,[n, 1]),2);

% Compute predicted covariance
P = Q;
for i = 1:N
P = P + W(i)*(fSP(:,i)-x)*(fSP(:,i)-x)';
end

% Make sure P is symmetric
P = 0.5*(P + P');

end

function [x, P] = updateMeanAndCovWithSigmaPoints(x, P, y, s, h, R, sigmaPoints, type)
%
%UPDATEGAUSSIANWITHSIGMAPOINTS computes the updated mean and covariance
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement
%   s           [2 x 1] sensor position
%   h           measurement model function handle
%   R           [m x m] measurement noise covariance matrix
%
%Output:
%   x           [n x 1] Updated mean
%   P           [n x n] Updated covariance
%

% Compute sigma points
[SP,W] = sigmaPoints(x, P, type);

% Dimension of measurement
m = size(R,1);

% Dimension of state and number of sigma points
[n, N] = size(SP);

% Predicted measurement
yhat = zeros(m,1);
hSP = zeros(m,N);
for i = 1:N
[hSP(:,i),~] = h(SP(:,i),s);
yhat = yhat + W(i)*hSP(:,i);
end

% Cross covariance and innovation covariance
Pxy = zeros(n,m);
S = R;
for i=1:N
Pxy = Pxy + W(i)*(SP(:,i)-x)*(hSP(:,i)-yhat)';
S = S + W(i)*(hSP(:,i)-yhat)*(hSP(:,i)-yhat)';
end

% Ensure symmetry
S = 0.5*(S+S');

% Updated mean
x = x+Pxy*(S\(y-yhat));
P = P - Pxy*(S\(Pxy'));

% Ensure symmetry
P = 0.5*(P+P');

end