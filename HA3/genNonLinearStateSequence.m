function X = genNonLinearStateSequence(x_0, P_0, f, Q, N)
  %GENNONLINEARSTATESEQUENCE generates an N+1-long sequence of states using a
  %    Gaussian prior and a nonlinear Gaussian process model
  %
  %Input:
  %   x_0         [n x 1] Prior mean
  %   P_0         [n x n] Prior covariance
  %   f           Motion model function handle
  %               [fx,Fx]=f(x)
  %               Takes as input x (state),
  %               Returns fx and Fx, motion model and Jacobian evaluated at x
  %               All other model parameters, such as sample time T,
  %               must be included in the function
  %   Q           [n x n] Process noise covariance
  %   N           [1 x 1] Number of states to generate
  %
  %Output:
  %   X           [n x N+1] State vector sequence
  %

  X = zeros(size(x_0,1), N+1);
  X(:,1) = mvnrnd(x_0, P_0);

  for i = 2:N+1
      X(:,i) = f(X(:,i-1))' + mvnrnd(zeros(size(Q,1),1), Q);
  end

end