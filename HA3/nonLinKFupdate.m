function [x, P] = nonLinKFupdate(x, P, y, h, R, type)
%NONLINKFUPDATE calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   y           [m x 1] measurement vector
%   h           Measurement model function handle
%               [hx,Hx]=h(x)
%               Takes as input x (state),
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%               Function must include all model parameters for the particular model,
%               such as sensor position for some models.
%   R           [m x m] Measurement noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] updated state mean
%   P           [n x n] updated state covariance
%

n = size(x,1);

    switch type
        case 'EKF'
            [x_mes, Jac] = h(x);

            S = Jac*P*Jac' + R;
            K = P*Jac'/S;

            x = x + K*(y-x_mes);
            P = P - K*S*K';

        case 'UKF'
            [Chis, Ws] = sigmaPoints(x, P, 'UKF');

            y_pred = 0;
            for i=1:2*n+1
                y_pred = y_pred + h(Chis(:,i))*Ws(i);
            end

            S = R;
            for i=1:2*n+1
                S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
            end

            Pxy = 0;
            for i=1:2*n+1
                Pxy = Pxy + (Chis(:,i) - x)*(h(Chis(:,i)) - y)'*Ws(i);
            end

            x = x + Pxy/S*(y-y_pred);
            P = P - Pxy/S*Pxy';

            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end

        case 'CKF'
            [Chis, Ws] = sigmaPoints(x, P, 'CKF');

            y_pred = 0;
            for i=1:2*n
                y_pred = y_pred + h(Chis(:,i))*Ws(i);
            end

            S = R;
            for i=1:2*n
                S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
            end

            Pxy = 0;
            for i=1:2*n
                Pxy = Pxy + (Chis(:,i) - x)*(h(Chis(:,i)) - y)'*Ws(i);
            end

            x = x + Pxy/S*(y-y_pred);
            P = P - Pxy/S*Pxy';

        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end

