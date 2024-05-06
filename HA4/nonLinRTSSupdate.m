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