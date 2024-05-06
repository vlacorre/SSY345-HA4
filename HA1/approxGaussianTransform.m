% Function retrieved from part HA1.1.3.
function [mu_y, Sigma_y, y_s] = approxGaussianTransform(mu_x, Sigma_x, f, N)
%approxGaussianTransform takes a Gaussian density and a transformation 
%function and calculates the mean and covariance of the transformed density.
%
%Inputs
%   MU_X        [m x 1] Expected value of x.
%   SIGMA_X     [m x m] Covariance of x.
%   F           [Function handle] Function which maps a [m x 1] dimensional
%               vector into another vector of size [n x 1].
%   N           Number of samples to draw. Default = 5000.
%
%Output
%   MU_Y        [n x 1] Approximated mean of y.
%   SIGMA_Y     [n x n] Approximated covariance of y.
%   ys          [n x N] Samples propagated through f

if nargin < 4
    N = 5000;
end

% Draw N random vectors with a normal distribution
x_s = mvnrnd(mu_x, Sigma_x, N)'; % This transpose results in size [m x N]

% x_s(:, i) is of size [m x 1]

% Get the size of the output vector of f to determine n
n = size(f(x_s(:, 1)), 1);

y_s = zeros(n, N);
for i = 1:N
    y_s(:, i) = f(x_s(:, i));
end

% Note: MATLAB default functions operate on rows. This is why we need to tweak the way we call mean() and cov().
mu_y = mean(y_s, 2); % Get the mean along the 2nd dimension of y_s

Sigma_y = cov(y_s');
end