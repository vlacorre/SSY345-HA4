function [x, P] = nonLinKFprediction(x, P, f, Q, type)
%NONLINKFPREDICTION calculates mean and covariance of predicted state
%   density using a non-linear Gaussian model.
%
%Input:
%   x           [n x 1] Prior mean
%   P           [n x n] Prior covariance
%   f           Motion model function handle
%               [fx,Fx]=f(x) 
%               Takes as input x (state), 
%               Returns fx and Fx, motion model and Jacobian evaluated at x
%               All other model parameters, such as sample time T,
%               must be included in the function
%   Q           [n x n] Process noise covariance
%   type        String that specifies the type of non-linear filter
%
%Output:
%   x           [n x 1] predicted state mean
%   P           [n x n] predicted state covariance
%

n = size(x,1);

    switch type
        case 'EKF'
            [x, Jac] = f(x);
            P = Jac*P*Jac' + Q;
            
        case 'UKF'
            [Chis, Ws] = sigmaPoints(x, P, 'UKF');

            x = 0;
            for i=1:2*n+1
                x = x + f(Chis(:,i))*Ws(i);
            end
            
            P = Q;
            for i=1:2*n+1
                P = P + (f(Chis(:,i)) - x)*(f(Chis(:,i)) - x)'*Ws(i);
            end
            
            % Make sure the covariance matrix is semi-definite
            if min(eig(P))<=0
                [v,e] = eig(P, 'vector');
                e(e<0) = 1e-4;
                P = v*diag(e)/v;
            end
                
        case 'CKF'
            [Chis, Ws] = sigmaPoints(x, P, 'CKF');

            x = 0;
            for i=1:2*n
                x = x + f(Chis(:,i))*Ws(i);
            end
            
            P = Q;
            for i=1:2*n
                P = P + (f(Chis(:,i)) - x)*(f(Chis(:,i)) - x)'*Ws(i);
            end
            
        otherwise
            error('Incorrect type of non-linear Kalman filter')
    end

end