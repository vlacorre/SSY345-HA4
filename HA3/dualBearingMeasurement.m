function [hx, Hx] = dualBearingMeasurement(x, s1, s2)
%DUALBEARINGMEASUREMENT calculates the bearings from two sensors, located in
%s1 and s2, to the position given by the state vector x. Also returns the
%Jacobian of the model at x.
%
%Input:
%   x           [n x 1] State vector, the two first element are 2D position
%   s1          [2 x 1] Sensor position (2D) for sensor 1
%   s2          [2 x 1] Sensor position (2D) for sensor 2
%
%Output:
%   hx          [2 x 1] measurement vector
%   Hx          [2 x n] measurement model Jacobian
%
% NOTE: the measurement model assumes that in the state vector x, the first
% two states are X-position and Y-position.

X = x(1);
Y = x(2);

s1_X = s1(1);
s1_Y = s1(2);

s2_X = s2(1);
s2_Y = s2(2);


hx_1 = atan2(Y - s1_Y, X - s1_X);
hx_2 = atan2(Y - s2_Y, X - s2_X);


hx = [hx_1; hx_2];

n = size(x,1);

Hx = zeros(2, n);
y1 = Y-s1_Y;
x1 = X-s1_X;
Hx(1,1) = -y1/(x1^2 + y1^2);
Hx(1,2) = x1/(x1^2 + y1^2);

y2 = Y-s2_Y;
x2 = X-s2_X;
Hx(2,1) = -y2/(x2^2 + y2^2);
Hx(2,2) = x2/(x2^2 + y2^2);

end