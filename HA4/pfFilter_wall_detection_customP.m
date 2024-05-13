function [xfp, Pfp, Xp, Wp] = pfFilter_wall_detection_customP(x_0, P_0, Y, proc_f, proc_Q, meas_h, meas_R, ...
    N, bResample, plotFunc)
    %PFFILTER Filters measurements Y using the SIS or SIR algorithms and a
    % state-space model.
    %
    % Input:
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   Y           [m x K] Measurement sequence to be filtered
    %   proc_f      Handle for process function f(x_k-1)
    %   proc_Q      [n x n] process noise covariance
    %   meas_h      Handle for measurement model function h(x_k)
    %   meas_R      [m x m] measurement noise covariance
    %   N           Number of particles
    %   bResample   boolean false - no resampling, true - resampling
    %   plotFunc    Handle for plot function that is called when a filter
    %               recursion has finished.
    % Output:
    %   xfp         [n x K] Posterior means of particle filter
    %   Pfp         [n x n x K] Posterior error covariances of particle filter
    %   Xp          [n x N x K] Non-resampled Particles for posterior state distribution in times 1:K
    %   Wp          [N x K] Non-resampled weights for posterior state x in times 1:K

    n = size(x_0, 1);
    K = size(Y, 2);

    xfp = zeros(n, K);
    Pfp = zeros(n, n, K);
    % Only store the particles if the function is called with more than 2 output arguments.
    if nargout > 2
        Xp = zeros(n, N, K);
        Wp = zeros(N, K);
    end

    X_k = zeros(n, N);
    W_k = zeros(1, N);
    for i = 1:N
        X_k(:,i) = mvnrnd(proc_f(x_0), P_0);
        W_k(i) = 1/N;
    end

    for k = 1:K
        X_kmin1 = X_k;
        [X_k, W_k] = pfFilterStep(X_k, W_k, Y(:,k), proc_f, proc_Q, meas_h, meas_R);
        if bResample
            [X_k, W_k, j] = resampl(X_k, W_k);
        else
            j = 1:N;
        end

        % [bool] = isOnRoad(x,y);
        % Resample particles that are not on road
        for i = 1:N
            original = X_k(:,i);
            count = 0;
            while ~isOnRoad(X_k(1,i), X_k(2,i)) && count < 100
                count = count + 1;
                X_k(:,i) = mvnrnd(proc_f(original), diag([5 5 0.1 0.1]));
            end
        end

        xfp(:,k) = X_k*W_k';
        Pfp(:,:,k) = (X_k-xfp(:,k)) * ((X_k-xfp(:,k))'.* W_k');

        if nargin(plotFunc) == 8
            ax = [-10 10 -10 10]; % ax = [xmin xmax ymin ymax];
            timeStepsToPlot = [2 15 29];
            if any(k == timeStepsToPlot)
                plotFunc(k, X_k, W_k, xfp, Pfp, bResample, 1, ax);
            end
        end

        if nargin(plotFunc) == 5
            plotFunc(k, X_k, X_kmin1, W_k, j);
        end

        % Only output the particles if the function is called with more than 2 output arguments.
        if nargout > 2
            Xp(:,:,k) = X_k;
            Wp(:,k) = W_k;
        end

    end

end

function [Xr, Wr, j] = resampl(X, W)
    %RESAMPLE Resample particles and output new particles and weights.
    % resampled particles.
    %
    %   if old particle vector is x, new particles x_new is computed as x(:,j)
    %
    % Input:
    %   X   [n x N] Particles, each column is a particle.
    %   W   [1 x N] Weights, corresponding to the samples
    %
    % Output:
    %   Xr  [n x N] Resampled particles, each corresponding to some particle
    %               from old weights.
    %   Wr  [1 x N] New weights for the resampled particles.
    %   j   [1 x N] vector of indices refering to vector of old particles
    Wc = cumsum(W);

    Xr = zeros(size(X));
    N = size(X,2);
    Wr = 1/N*ones(1,N);

    j = zeros(1,N);
    for i = 1:N
        u = rand;
        % Select the particle that corresponds to this weight zone
        for k = 1:N
            if u < Wc(k)
                Xr(:,i) = X(:,k);
                j(i) = k;
                break;
            end
        end
    end
end


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