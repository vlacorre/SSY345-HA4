function [SP,W] = sigmaPoints(x, P, type)
  % SIGMAPOINTS computes sigma points, either using unscented transform or
  % using cubature.
  %
  %Input:
  %   x           [n x 1] Prior mean
  %   P           [n x n] Prior covariance
  %
  %Output:
  %   SP          [n x 2n+1] UKF, [n x 2n] CKF. Matrix with sigma points
  %   W           [1 x 2n+1] UKF, [1 x 2n] UKF. Vector with sigma point weights
  %
  n = size(x,1);
  P_sqrt = sqrtm(P);

  switch type
    case 'UKF'
      W0 = 1 - n/3;
      P_sqrt = sqrtm(P);

      SP = zeros(n, 2*n+1);
      W = zeros(1, 2*n+1);

      SP(:, 1) = x;
      W(1) = W0;

      for i = 1:n
          SP(:, 1+i) = x + sqrt(n/(1-W0))*P_sqrt(:,i);
          W(1+i) = (1 - W0)/(2*n);
      end
      for i = 1:n
          SP(:, 1+n+i) = x - sqrt(n/(1-W0))*P_sqrt(:,i);
          W(1+n+i) = (1 - W0)/(2*n);
      end

    case 'CKF'
      SP = zeros(n, 2*n);
      W = zeros(1, 2*n);

      for i = 1:n
          SP(:, i) = x + sqrt(n)*P_sqrt(:,i);
          W(i) = 1/(2*n);
      end

      for i = 1:n
          SP(:, n+i) = x - sqrt(n)*P_sqrt(:,i);
          W(n+i) = 1/(2*n);
      end

    otherwise
      error('Incorrect type of sigma point')
  end

end