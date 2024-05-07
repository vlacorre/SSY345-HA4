
function [X_k, W_k] = pfFilterStep(X_kmin1, W_kmin1, yk, proc_f, proc_Q, meas_h, meas_R)
    %PFFILTERSTEP Compute one filter step of a SIS/SIR particle filter.
    %
    % Input:
    %   X_kmin1     [n x N] Particles for state x in time k-1
    %   W_kmin1     [1 x N] Weights for state x in time k-1
    %   y_k         [m x 1] Measurement vector for time k
    %   proc_f      Handle for process function f(x_k-1)
    %   proc_Q      [n x n] process noise covariance
    %   meas_h      Handle for measurement model function h(x_k)
    %   meas_R      [m x m] measurement noise covariance
    %
    % Output:
    %   X_k         [n x N] Particles for state x in time k
    %   W_k         [1 x N] Weights for state x in time k

    N = size(X_kmin1, 2);
    X_k = zeros(size(X_kmin1));
    W_k = zeros(1, N);

    for i = 1:N
        % Proposal density is the process model
        X_k(:,i) = mvnrnd(proc_f(X_kmin1(:,i)), proc_Q);
        % Likelihood is the measurement model given the particle
        W_k(i) = W_kmin1(i) * mvnpdf(yk, meas_h(X_k(:,i)), meas_R);
    end

    % Normalize weights
    sum_Wk = sum(W_k);
    W_k = W_k / sum_Wk;

end