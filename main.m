set(0, 'DefaultAxesXGrid', 'on');  % Turn on grid for x-axis
set(0, 'DefaultAxesYGrid', 'on');  % Turn on grid for y-axis

set(0, 'DefaultLegendFontSizeMode', 'manual');
set(0, 'DefaultLegendFontSize', 14);
set(0, 'DefaultLegendLocation', 'best');


set(0, 'DefaultAxesFontSize', 14);  % Increase the default font size for axes
set(0, 'DefaultTextFontSize', 14);  % Increase the default font size for text


set(0, 'DefaultAxesFontSize', 30);  % Increase the default font size for axes
set(0, 'DefaultLineLineWidth', 3);  % Increase the default line width
set(0, 'DefaultTextFontSize', 26);  % Increase the default font size for text
set(0, 'DefaultAxesTitleFontSizeMultiplier', 1);  % Increase the default font size for axes titles
set(0, 'DefaultAxesLabelFontSizeMultiplier', 1.1);  % Increase the default font size for axes labels
set(0, 'DefaultFigurePosition', [0, -1500, 1900, 1500]); % [x, y, width, height]
set(0, 'DefaultLineMarkersize', 10);

% Increase legend box
set(0, 'DefaultLegendBox', 'on');
set(0, 'DefaultLegendFontSizeMode', 'manual');
set(0, 'DefaultLegendFontSize', 20);
set(0, 'DefaultLegendLocation', 'best');

clc;  clear all;  close all;

% Delete old pictures
delete('Images/*.eps');

addpath('HA1');
addpath('HA2');
addpath('HA3');

rng(1); % Set random seed

scenario1 = true;
scenario2 = true;
scenario3 = true;

if scenario1
%% Scenario 1: Approximations of mean and covariance
%% a) Measurement samples

% Different state densities
x_1_hat = [125; 125];
x_2_hat = [-25; 125];
x_3_hat = [60; 60];
P123 = diag([10, 5].^2);

s1 = [0; 100];
s2 = [100; 0];
h = @(x) dualBearingMeasurement(x, s1, s2); % Measurement model
sigma_phi = 0.1*pi/180; % Standard deviation of the measurement noise
R = diag(([1, 1].*sigma_phi).^2); % Measurement noise covariance

% Approximation of the mean and covariance of the measurement sequence for different state densities
% Density 1
figure;
N = 10000;
X_true = mvnrnd(x_1_hat, P123, N)'; % Generate the true state sequence
Y = genNonLinearMeasurementSequence(X_true, h, R); % Generate the measurement sequence
y_bar = mean(Y, 2);
var_y_bar = cov(Y');

plot(Y(1,:), Y(2,:), '.', 'Color', [.4 .4 .4]);
hold on;
plot(y_bar(1), y_bar(2), 'rx');
xy_1 = sigmaEllipse2D(y_bar, var_y_bar, 1, N);
xy_3 = sigmaEllipse2D(y_bar, var_y_bar, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'r-');
plot(xy_3(1,:), xy_3(2,:), 'r--');
hold off;
xlabel('y_1');
ylabel('y_2');
legend('Measurement samples', 'Mean', '1\sigma level', '3\sigma level', 'Location', 'eastoutside');

print('Images/1_a_approx_density1.eps', '-depsc');

% Density 2
figure;
N = 10000;
X_true = mvnrnd(x_2_hat, P123, N)'; % Generate the true state sequence
Y = genNonLinearMeasurementSequence(X_true, h, R); % Generate the measurement sequence
y_bar = mean(Y, 2);
var_y_bar = cov(Y');

plot(Y(1,:), Y(2,:), '.', 'Color', [.4 .4 .4]);
hold on;
plot(y_bar(1), y_bar(2), 'rx');
xy_1 = sigmaEllipse2D(y_bar, var_y_bar, 1, N);
xy_3 = sigmaEllipse2D(y_bar, var_y_bar, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'r-');
plot(xy_3(1,:), xy_3(2,:), 'r--');
hold off;
xlabel('y_1');
ylabel('y_2');
legend('Measurement samples', 'Mean', '1\sigma level', '3\sigma level', 'Location', 'eastoutside');

print('Images/1_a_approx_density2.eps', '-depsc');

% Density 3
figure;
N = 10000;
X_true = mvnrnd(x_3_hat, P123, N)'; % Generate the true state sequence
Y = genNonLinearMeasurementSequence(X_true, h, R); % Generate the measurement sequence
y_bar = mean(Y, 2);
var_y_bar = cov(Y');

plot(Y(1,:), Y(2,:), '.', 'Color', [.4 .4 .4]);
hold on;
plot(y_bar(1), y_bar(2), 'rx');
xy_1 = sigmaEllipse2D(y_bar, var_y_bar, 1, N);
xy_3 = sigmaEllipse2D(y_bar, var_y_bar, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'r-');
plot(xy_3(1,:), xy_3(2,:), 'r--');
hold off;
xlabel('y_1');
ylabel('y_2');
legend('Measurement samples', 'Mean', '1\sigma level', '3\sigma level', 'Location', 'eastoutside');

print('Images/1_a_approx_density3.eps', '-depsc');

%% b) Analytical approximations
% Density 1
% EKF
[hx, Hx] = h(x_1_hat);
y_1_bar_ekf = hx
var_y_1_hat_ekf = Hx*P123*Hx' + R

% UKF
n = size(x_1_hat,1);
[Chis, Ws] = sigmaPoints(x_1_hat, P123, 'UKF');
y_pred = 0;
for i=1:2*n+1
    y_pred = y_pred + h(Chis(:,i))*Ws(i);
end
S = R;
for i=1:2*n+1
    S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
end
y_1_bar_ukf = y_pred
var_y_1_hat_ukf = S

% CKF
n = size(x_1_hat,1);
[Chis, Ws] = sigmaPoints(x_1_hat, P123, 'CKF');
y_pred = 0;
for i=1:2*n
    y_pred = y_pred + h(Chis(:,i))*Ws(i);
end
S = R;
for i=1:2*n
    S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
end
y_1_bar_ckf = y_pred
var_y_1_hat_ckf = S

% Density 2
% EKF
[hx, Hx] = h(x_2_hat);
y_2_bar_ekf = hx
var_y_2_hat_ekf = Hx*P123*Hx' + R

% UKF
n = size(x_2_hat,1);
[Chis, Ws] = sigmaPoints(x_2_hat, P123, 'UKF');
y_pred = 0;
for i=1:2*n+1
    y_pred = y_pred + h(Chis(:,i))*Ws(i);
end
S = R;
for i=1:2*n+1
    S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
end
y_2_bar_ukf = y_pred
var_y_2_hat_ukf = S

% CKF
n = size(x_2_hat,1);
[Chis, Ws] = sigmaPoints(x_2_hat, P123, 'CKF');
y_pred = 0;
for i=1:2*n
    y_pred = y_pred + h(Chis(:,i))*Ws(i);
end
S = R;
for i=1:2*n
    S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
end
y_2_bar_ckf = y_pred
var_y_2_hat_ckf = S

% Density 3
% EKF
[hx, Hx] = h(x_3_hat);
y_3_bar_ekf = hx
var_y_3_hat_ekf = Hx*P123*Hx' + R

% UKF
n = size(x_3_hat,1);
[Chis, Ws] = sigmaPoints(x_3_hat, P123, 'UKF');
y_pred = 0;
for i=1:2*n+1
    y_pred = y_pred + h(Chis(:,i))*Ws(i);
end
S = R;
for i=1:2*n+1
    S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
end
y_3_bar_ukf = y_pred
var_y_3_hat_ukf = S

% CKF
n = size(x_3_hat,1);
[Chis, Ws] = sigmaPoints(x_3_hat, P123, 'CKF');
y_pred = 0;
for i=1:2*n
    y_pred = y_pred + h(Chis(:,i))*Ws(i);
end
S = R;
for i=1:2*n
    S = S + (h(Chis(:,i)) - y_pred)*(h(Chis(:,i)) - y_pred)'*Ws(i);
end
y_3_bar_ckf = y_pred
var_y_3_hat_ckf = S


%% c) Comparison

% Density 1
figure;
N = 10000;
X_true = mvnrnd(x_1_hat, P123, N)'; % Generate the true state sequence
Y = genNonLinearMeasurementSequence(X_true, h, R); % Generate the measurement sequence
y_bar = mean(Y, 2);
var_y_bar = cov(Y');

plot(Y(1,:), Y(2,:), '.', 'Color', [.4 .4 .4]);
hold on;
plot(y_bar(1), y_bar(2), 'rx');
xy_1 = sigmaEllipse2D(y_bar, var_y_bar, 1, N);
xy_3 = sigmaEllipse2D(y_bar, var_y_bar, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'r-');
plot(xy_3(1,:), xy_3(2,:), 'r--');

plot(y_1_bar_ekf(1), y_1_bar_ekf(2), 'bx');
xy_1 = sigmaEllipse2D(y_1_bar_ekf, var_y_1_hat_ekf, 1, N);
xy_3 = sigmaEllipse2D(y_1_bar_ekf, var_y_1_hat_ekf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'b-');
plot(xy_3(1,:), xy_3(2,:), 'b--');


plot(y_1_bar_ukf(1), y_1_bar_ukf(2), 'mx');
xy_1 = sigmaEllipse2D(y_1_bar_ukf, var_y_1_hat_ukf, 1, N);
xy_3 = sigmaEllipse2D(y_1_bar_ukf, var_y_1_hat_ukf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'm-');
plot(xy_3(1,:), xy_3(2,:), 'm--');


plot(y_1_bar_ckf(1), y_1_bar_ckf(2), 'gx');
xy_1 = sigmaEllipse2D(y_1_bar_ckf, var_y_1_hat_ckf, 1, N);
xy_3 = sigmaEllipse2D(y_1_bar_ckf, var_y_1_hat_ckf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'g-');
plot(xy_3(1,:), xy_3(2,:), 'g--');
hold off;
xlabel('y_1');
ylabel('y_2');
legend('Measurement samples', 'Mean', '1\sigma level', '3\sigma level', 'EKF', 'EKF 1\sigma level', 'EKF 3\sigma level', 'UKF', 'UKF 1\sigma level', 'UKF 3\sigma level', 'CKF', 'CKF 1\sigma level', 'CKF 3\sigma level');

print('Images/1_c_analytic1.eps', '-depsc');

% Density 2
figure;
N = 10000;
X_true = mvnrnd(x_2_hat, P123, N)'; % Generate the true state sequence
Y = genNonLinearMeasurementSequence(X_true, h, R); % Generate the measurement sequence
y_bar = mean(Y, 2);
var_y_bar = cov(Y');

plot(Y(1,:), Y(2,:), '.', 'Color', [.4 .4 .4]);
hold on;
plot(y_bar(1), y_bar(2), 'rx');
xy_1 = sigmaEllipse2D(y_bar, var_y_bar, 1, N);
xy_3 = sigmaEllipse2D(y_bar, var_y_bar, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'r-');
plot(xy_3(1,:), xy_3(2,:), 'r--');

plot(y_2_bar_ekf(1), y_2_bar_ekf(2), 'bx');
xy_1 = sigmaEllipse2D(y_2_bar_ekf, var_y_2_hat_ekf, 1, N);
xy_3 = sigmaEllipse2D(y_2_bar_ekf, var_y_2_hat_ekf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'b-');
plot(xy_3(1,:), xy_3(2,:), 'b--');


plot(y_2_bar_ukf(1), y_2_bar_ukf(2), 'mx');
xy_1 = sigmaEllipse2D(y_2_bar_ukf, var_y_2_hat_ukf, 1, N);
xy_3 = sigmaEllipse2D(y_2_bar_ukf, var_y_2_hat_ukf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'm-');
plot(xy_3(1,:), xy_3(2,:), 'm--');


plot(y_2_bar_ckf(1), y_2_bar_ckf(2), 'gx');
xy_1 = sigmaEllipse2D(y_2_bar_ckf, var_y_2_hat_ckf, 1, N);
xy_3 = sigmaEllipse2D(y_2_bar_ckf, var_y_2_hat_ckf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'g-');
plot(xy_3(1,:), xy_3(2,:), 'g--');
hold off;
xlabel('y_1');
ylabel('y_2');
legend('Measurement samples', 'Mean', '1\sigma level', '3\sigma level', 'EKF', 'EKF 1\sigma level', 'EKF 3\sigma level', 'UKF', 'UKF 1\sigma level', 'UKF 3\sigma level', 'CKF', 'CKF 1\sigma level', 'CKF 3\sigma level');

print('Images/1_c_analytic2.eps', '-depsc');

% Density 3
figure;
N = 10000;
X_true = mvnrnd(x_3_hat, P123, N)'; % Generate the true state sequence
Y = genNonLinearMeasurementSequence(X_true, h, R); % Generate the measurement sequence
y_bar = mean(Y, 2);
var_y_bar = cov(Y');

plot(Y(1,:), Y(2,:), '.', 'Color', [.4 .4 .4]);
hold on;
plot(y_bar(1), y_bar(2), 'rx');
xy_1 = sigmaEllipse2D(y_bar, var_y_bar, 1, N);
xy_3 = sigmaEllipse2D(y_bar, var_y_bar, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'r-');
plot(xy_3(1,:), xy_3(2,:), 'r--');

plot(y_3_bar_ekf(1), y_3_bar_ekf(2), 'bx');
xy_1 = sigmaEllipse2D(y_3_bar_ekf, var_y_3_hat_ekf, 1, N);
xy_3 = sigmaEllipse2D(y_3_bar_ekf, var_y_3_hat_ekf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'b-');
plot(xy_3(1,:), xy_3(2,:), 'b--');


plot(y_3_bar_ukf(1), y_3_bar_ukf(2), 'mx');
xy_1 = sigmaEllipse2D(y_3_bar_ukf, var_y_3_hat_ukf, 1, N);
xy_3 = sigmaEllipse2D(y_3_bar_ukf, var_y_3_hat_ukf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'm-');
plot(xy_3(1,:), xy_3(2,:), 'm--');


plot(y_3_bar_ckf(1), y_3_bar_ckf(2), 'gx');
xy_1 = sigmaEllipse2D(y_3_bar_ckf, var_y_3_hat_ckf, 1, N);
xy_3 = sigmaEllipse2D(y_3_bar_ckf, var_y_3_hat_ckf, 3, N);
plot(xy_1(1,:), xy_1(2,:), 'g-');
plot(xy_3(1,:), xy_3(2,:), 'g--');
hold off;
xlabel('y_1');
ylabel('y_2');
legend('Measurement samples', 'Mean', '1\sigma level', '3\sigma level', 'EKF', 'EKF 1\sigma level', 'EKF 3\sigma level', 'UKF', 'UKF 1\sigma level', 'UKF 3\sigma level', 'CKF', 'CKF 1\sigma level', 'CKF 3\sigma level');

print('Images/1_c_analytic3.eps', '-depsc');

end % scenario1

if scenario2

%% Scenario 2: Non-linear Kalman filtering
% Prior
x_0 = [0; 0; 20; 0; 5*pi/180];
P_0 = diag([10, 10, 2, pi/180, pi/180].^2);

% Motion model
T = 1; % Sampling time
f = @(x) coordinatedTurnMotion(x, T); % Motion model
Q123 = diag([0, 0, 1, 0, pi/180].^2); % Process noise covariance (Case 1, 2 and 3)


% Sensor
s1 = [-200; 100];
s2 = [-200; -100];
h = @(x) dualBearingMeasurement(x, s1, s2); % Measurement model
R1 = diag([2*pi/180, 2*pi/180].^2); % Measurement noise covariance (Case 1)
R2 = diag([2*pi/180, 0.1*pi/180].^2); % Measurement noise covariance (Case 2)
R3 = diag([0.1*pi/180, 0.1*pi/180].^2); % Measurement noise covariance (Case 3)

%% a) Generate true state and measurement sequences (Case 1)
N = 100;
X_true = genNonLinearStateSequence(x_0, P_0, f, Q123, N);

Y = genNonLinearMeasurementSequence(X_true, h, R1);

sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for i=1:N
    [sensorPositions_x(i), sensorPositions_y(i)] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
end
[xf_ekf, Pf_ekf, xp_ekf, Pp_ekf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R1, 'EKF');
[xf_ukf, Pf_ukf, xp_ukf, Pp_ukf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R1, 'UKF');
[xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R1, 'CKF');



% Plot the true and estimated states
figure;
plot(X_true(1,:), X_true(2,:), 'r--x');
hold on;
plot(xf_ekf(1,:), xf_ekf(2,:), 'b--x');
plot(xf_ukf(1,:), xf_ukf(2,:), 'm--x');
plot(xf_ckf(1,:), xf_ckf(2,:), 'g--x');
plot(X_true(1,1), X_true(2,1), 'ro'); % Add a circle at the first point of X_true
plot(xf_ekf(1,1), xf_ekf(2,1), 'black o'); % Add a circle at the first point of X_true
plot(sensorPositions_x, sensorPositions_y, 'black x');
hold off;
xlabel('x');
ylabel('y');
legend('True state', 'EKF estimates', 'UKF estimates', 'CKF estimates', 'True initial position', 'Predicted initial position', 'Sensor positions', 'Location', 'eastoutside');
xlim([-1200, 260]);
ylim([-100, 1400]);

print('Images/2_a_case1_allTypes.eps', '-depsc');


% Plot the true and estimated states
figure;
plot(X_true(1,:), X_true(2,:), 'r--x');
hold on;
plot(xf_ekf(1,:), xf_ekf(2,:), 'b--x');
plot(xf_ukf(1,:), xf_ukf(2,:), 'm--x');
plot(xf_ckf(1,:), xf_ckf(2,:), 'g--x');
plot(X_true(1,1), X_true(2,1), 'ro'); % Add a circle at the first point of X_true
plot(xf_ekf(1,1), xf_ekf(2,1), 'black o'); % Add a circle at the first point of X_true
plot(sensorPositions_x, sensorPositions_y, 'black x');

% Plot the 3-sigma level curve around every 5th point
for i=5:5:N
    xy_3_ekf = sigmaEllipse2D(xf_ekf(1:2,i), Pf_ekf(1:2,1:2,i), 3, 100);
    xy_3_ukf = sigmaEllipse2D(xf_ukf(1:2,i), Pf_ukf(1:2,1:2,i), 3, 100);
    xy_3_ckf = sigmaEllipse2D(xf_ckf(1:2,i), Pf_ckf(1:2,1:2,i), 3, 100);
    plot(xy_3_ekf(1,:), xy_3_ekf(2,:), 'r-');
    plot(xy_3_ukf(1,:), xy_3_ukf(2,:), 'm-');
    plot(xy_3_ckf(1,:), xy_3_ckf(2,:), 'g-');
end

hold off;
xlabel('x');
ylabel('y');
legend('True state', 'EKF estimates', 'UKF estimates', 'CKF estimates', 'True initial position', 'Predicted initial position', 'Sensor positions', 'EKF Uncertainty', 'UKF Uncertainty', 'CKF Uncertainty', 'Location', 'eastoutside');
xlim([-1200, 260]);
ylim([-100, 1400]);

print('Images/2_a_case1_allTypes_uncertainties.eps', '-depsc');

%% b) Case 2 and 3
% Case 2
Y = genNonLinearMeasurementSequence(X_true, h, R2);

sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for i=1:N
    [sensorPositions_x(i), sensorPositions_y(i)] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
end
[xf_ekf, Pf_ekf, xp_ekf, Pp_ekf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R2, 'EKF');
[xf_ukf, Pf_ukf, xp_ukf, Pp_ukf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R2, 'UKF');
[xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R2, 'CKF');

% Plot the true and estimated states
figure;
plot(X_true(1,:), X_true(2,:), 'r--x');
hold on;
plot(xf_ekf(1,:), xf_ekf(2,:), 'b--x');
plot(xf_ukf(1,:), xf_ukf(2,:), 'm--x');
plot(xf_ckf(1,:), xf_ckf(2,:), 'g--x');
plot(X_true(1,1), X_true(2,1), 'ro'); % Add a circle at the first point of X_true
plot(xf_ekf(1,1), xf_ekf(2,1), 'black o'); % Add a circle at the first point of X_true
plot(sensorPositions_x, sensorPositions_y, 'black x');

% Plot the 3-sigma level curve around every 5th point
for i=5:5:N
    xy_3_ekf = sigmaEllipse2D(xf_ekf(1:2,i), Pf_ekf(1:2,1:2,i), 3, 100);
    xy_3_ukf = sigmaEllipse2D(xf_ukf(1:2,i), Pf_ukf(1:2,1:2,i), 3, 100);
    xy_3_ckf = sigmaEllipse2D(xf_ckf(1:2,i), Pf_ckf(1:2,1:2,i), 3, 100);
    plot(xy_3_ekf(1,:), xy_3_ekf(2,:), 'r-');
    plot(xy_3_ukf(1,:), xy_3_ukf(2,:), 'm-');
    plot(xy_3_ckf(1,:), xy_3_ckf(2,:), 'g-');
end
hold off;
xlabel('x');
ylabel('y');
legend('True state', 'EKF estimates', 'UKF estimates', 'CKF estimates', 'True initial position', 'Predicted initial position', 'Sensor positions', 'EKF Uncertainty', 'UKF Uncertainty', 'CKF Uncertainty', 'Location', 'eastoutside');
xlim([-1200, 260]);
ylim([-100, 1400]);

print('Images/2_b_case2_allTypes_uncertainties.eps', '-depsc');

% Case 3
Y = genNonLinearMeasurementSequence(X_true, h, R3);
sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for i=1:N
    [sensorPositions_x(i), sensorPositions_y(i)] = getPosFromMeasurement(Y(1,i), Y(2,i), s1, s2);
end
[xf_ekf, Pf_ekf, xp_ekf, Pp_ekf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R3, 'EKF');
[xf_ukf, Pf_ukf, xp_ukf, Pp_ukf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R3, 'UKF');
[xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y, x_0, P_0, f, Q123, h, R3, 'CKF');
% Plot the true and estimated states
figure;
plot(X_true(1,:), X_true(2,:), 'r--x');
hold on;
plot(xf_ekf(1,:), xf_ekf(2,:), 'b--x');
plot(xf_ukf(1,:), xf_ukf(2,:), 'm--x');
plot(xf_ckf(1,:), xf_ckf(2,:), 'g--x');
plot(X_true(1,1), X_true(2,1), 'ro'); % Add a circle at the first point of X_true
plot(xf_ekf(1,1), xf_ekf(2,1), 'black o'); % Add a circle at the first point of X_true
plot(sensorPositions_x, sensorPositions_y, 'black x');

% Plot the 3-sigma level curve around every 5th point
for i=5:5:N
    xy_3_ekf = sigmaEllipse2D(xf_ekf(1:2,i), Pf_ekf(1:2,1:2,i), 3, 100);
    xy_3_ukf = sigmaEllipse2D(xf_ukf(1:2,i), Pf_ukf(1:2,1:2,i), 3, 100);
    xy_3_ckf = sigmaEllipse2D(xf_ckf(1:2,i), Pf_ckf(1:2,1:2,i), 3, 100);
    plot(xy_3_ekf(1,:), xy_3_ekf(2,:), 'r-');
    plot(xy_3_ukf(1,:), xy_3_ukf(2,:), 'm-');
    plot(xy_3_ckf(1,:), xy_3_ckf(2,:), 'g-');
end
hold off;
xlabel('x');
ylabel('y');
legend('True state', 'EKF estimates', 'UKF estimates', 'CKF estimates', 'True initial position', 'Predicted initial position', 'Sensor positions', 'EKF Uncertainty', 'UKF Uncertainty', 'CKF Uncertainty', 'Location', 'eastoutside');
xlim([-1200, 260]);
ylim([-100, 1400]);

print('Images/2_b_case3_allTypes_uncertainties.eps', '-depsc');

%% c) True vs EKF vs CKF
x_error_ekf_c1 = [];
x_error_ckf_c1 = [];
y_error_ekf_c1 = [];
y_error_ckf_c1 = [];

x_error_ekf_c2 = [];
x_error_ckf_c2 = [];
y_error_ekf_c2 = [];
y_error_ckf_c2 = [];

x_error_ekf_c3 = [];
x_error_ckf_c3 = [];
y_error_ekf_c3 = [];
y_error_ckf_c3 = [];

divergence_threshold = 100;

for i=1:150
    X_true = genNonLinearStateSequence(x_0, P_0, f, Q123, N);

    Y1 = genNonLinearMeasurementSequence(X_true, h, R1);
    [xf_ekf, Pf_ekf, xp_ekf, Pp_ekf] = nonLinearKalmanFilter(Y1, x_0, P_0, f, Q123, h, R1, 'EKF');
    [xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y1, x_0, P_0, f, Q123, h, R1, 'CKF');
    % Calculate estimation errors
    for j=1:N
        if abs(xf_ekf(1,j) - X_true(1,j)) <= divergence_threshold
            x_error_ekf_c1 = [x_error_ekf_c1, xf_ekf(1,j) - X_true(1,j)];
        end
        if abs(xf_ckf(1,j) - X_true(1,j)) <= divergence_threshold
            x_error_ckf_c1 = [x_error_ckf_c1, xf_ckf(1,j) - X_true(1,j)];
        end
        if abs(xf_ekf(2,j) - X_true(2,j)) <= divergence_threshold
            y_error_ekf_c1 = [y_error_ekf_c1, xf_ekf(2,j) - X_true(2,j)];
        end
        if abs(xf_ckf(2,j) - X_true(2,j)) <= divergence_threshold
            y_error_ckf_c1 = [y_error_ckf_c1, xf_ckf(2,j) - X_true(2,j)];
        end
    end
    
    Y2 = genNonLinearMeasurementSequence(X_true, h, R2);
    [xf_ekf, Pf_ekf, xp_ekf, Pp_ekf] = nonLinearKalmanFilter(Y2, x_0, P_0, f, Q123, h, R2, 'EKF');
    [xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y2, x_0, P_0, f, Q123, h, R2, 'CKF');
    % Calculate estimation errors
    for j=1:N
        if abs(xf_ekf(1,j) - X_true(1,j)) <= divergence_threshold
            x_error_ekf_c2 = [x_error_ekf_c2, xf_ekf(1,j) - X_true(1,j)];
        end
        if abs(xf_ckf(1,j) - X_true(1,j)) <= divergence_threshold
            x_error_ckf_c2 = [x_error_ckf_c2, xf_ckf(1,j) - X_true(1,j)];
        end
        if abs(xf_ekf(2,j) - X_true(2,j)) <= divergence_threshold
            y_error_ekf_c2 = [y_error_ekf_c2, xf_ekf(2,j) - X_true(2,j)];
        end
        if abs(xf_ckf(2,j) - X_true(2,j)) <= divergence_threshold
            y_error_ckf_c2 = [y_error_ckf_c2, xf_ckf(2,j) - X_true(2,j)];
        end
    end

    Y3 = genNonLinearMeasurementSequence(X_true, h, R3);
    [xf_ekf, Pf_ekf, xp_ekf, Pp_ekf] = nonLinearKalmanFilter(Y3, x_0, P_0, f, Q123, h, R3, 'EKF');
    [xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y3, x_0, P_0, f, Q123, h, R3, 'CKF');
    % Calculate estimation errors
    for j=1:N
        if abs(xf_ekf(1,j) - X_true(1,j)) <= divergence_threshold
            x_error_ekf_c3 = [x_error_ekf_c3, xf_ekf(1,j) - X_true(1,j)];
        end
        if abs(xf_ckf(1,j) - X_true(1,j)) <= divergence_threshold
            x_error_ckf_c3 = [x_error_ckf_c3, xf_ckf(1,j) - X_true(1,j)];
        end
        if abs(xf_ekf(2,j) - X_true(2,j)) <= divergence_threshold
            y_error_ekf_c3 = [y_error_ekf_c3, xf_ekf(2,j) - X_true(2,j)];
        end
        if abs(xf_ckf(2,j) - X_true(2,j)) <= divergence_threshold
            y_error_ckf_c3 = [y_error_ckf_c3, xf_ckf(2,j) - X_true(2,j)];
        end
    end

    
end

% Plot histograms (case 1)
figure;
subplot(2,1,1);
xlim([-divergence_threshold divergence_threshold]);
ylim([0 0.02]);
hold on;
histogram(x_error_ekf_c1, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(x_error_ekf_c1);
sigma = std(x_error_ekf_c1);
x = linspace(min(x_error_ekf_c1), max(x_error_ekf_c1), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);

histogram(x_error_ckf_c1, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(x_error_ckf_c1);
sigma = std(x_error_ckf_c1);
x = linspace(min(x_error_ckf_c1), max(x_error_ckf_c1), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'g-', 'LineWidth', 2);


hold off;
xlabel('Estimation Error (x)');
ylabel('Probability');
legend('EKF', 'EKF Gaussian Fit', 'CKF', 'CKF Gaussian Fit');

subplot(2,1,2);
xlim([-divergence_threshold divergence_threshold]);
ylim([0 0.02]);
hold on;
histogram(y_error_ekf_c1, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(y_error_ekf_c1);
sigma = std(y_error_ekf_c1);
x = linspace(min(y_error_ekf_c1), max(y_error_ekf_c1), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);

histogram(y_error_ckf_c1, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(y_error_ckf_c1);
sigma = std(y_error_ckf_c1);
x = linspace(min(y_error_ckf_c1), max(y_error_ckf_c1), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'g-', 'LineWidth', 2);

hold off;
xlabel('Estimation Error (y)');
ylabel('Probability');
legend('EKF', 'EKF Gaussian Fit', 'CKF', 'CKF Gaussian Fit');
print('Images/2_c_histograms_case1.eps', '-depsc');


% Plot histograms (case 2)
figure;
subplot(2,1,1);
xlim([-divergence_threshold divergence_threshold]);
ylim([0 0.02]);
hold on;
histogram(x_error_ekf_c2, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(x_error_ekf_c2);
sigma = std(x_error_ekf_c2);
x = linspace(min(x_error_ekf_c2), max(x_error_ekf_c2), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);

histogram(x_error_ckf_c2, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(x_error_ckf_c2);
sigma = std(x_error_ckf_c2);
x = linspace(min(x_error_ckf_c2), max(x_error_ckf_c2), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'g-', 'LineWidth', 2);

hold off;
xlabel('Estimation Error (x)');
ylabel('Probability');
legend('EKF', 'EKF Gaussian Fit', 'CKF', 'CKF Gaussian Fit');

subplot(2,1,2);
xlim([-divergence_threshold divergence_threshold]);
ylim([0 0.02]);
hold on;
histogram(y_error_ekf_c2, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(y_error_ekf_c2);
sigma = std(y_error_ekf_c2);
x = linspace(min(y_error_ekf_c2), max(y_error_ekf_c2), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);

histogram(y_error_ckf_c2, 'Normalization', 'pdf', 'NumBins', 100);
mu = mean(y_error_ckf_c2);
sigma = std(y_error_ckf_c2);
x = linspace(min(y_error_ckf_c2), max(y_error_ckf_c2), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'g-', 'LineWidth', 2);

hold off;
xlabel('Estimation Error (y)');
ylabel('Probability');
legend('EKF', 'EKF Gaussian Fit', 'CKF', 'CKF Gaussian Fit');
print('Images/2_c_histograms_case2.eps', '-depsc');


% Plot histograms (case 3)
figure;
subplot(2,1,1);
xlim([-divergence_threshold divergence_threshold]);
ylim([0 0.03]);
hold on;
histogram(x_error_ekf_c3, 'Normalization', 'pdf', 'NumBins', 200);
mu = mean(x_error_ekf_c3);
sigma = std(x_error_ekf_c3);
x = linspace(min(x_error_ekf_c3), max(x_error_ekf_c3), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);

histogram(x_error_ckf_c3, 'Normalization', 'pdf', 'NumBins', 200);
mu = mean(x_error_ckf_c3);
sigma = std(x_error_ckf_c3);
x = linspace(min(x_error_ckf_c3), max(x_error_ckf_c3), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'g-', 'LineWidth', 2);

hold off;
xlabel('Estimation Error (x)');
ylabel('Probability');
legend('EKF', 'EKF Gaussian Fit', 'CKF', 'CKF Gaussian Fit');

subplot(2,1,2);
xlim([-divergence_threshold divergence_threshold]);
ylim([0 0.03]);
hold on;
histogram(y_error_ekf_c3, 'Normalization', 'pdf', 'NumBins', 300);
mu = mean(y_error_ekf_c3);
sigma = std(y_error_ekf_c3);
x = linspace(min(y_error_ekf_c3), max(y_error_ekf_c3), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);
histogram(y_error_ckf_c3, 'Normalization', 'pdf', 'NumBins', 300);
mu = mean(y_error_ckf_c3);
sigma = std(y_error_ckf_c3);
x = linspace(min(y_error_ckf_c3), max(y_error_ckf_c3), 100);
y = normpdf(x, mu, sigma);
plot(x, y, 'g-', 'LineWidth', 2);
hold off;
xlabel('Estimation Error (y)');
ylabel('Probability');
legend('EKF', 'EKF Gaussian Fit', 'CKF', 'CKF Gaussian Fit');
print('Images/2_c_histograms_case3.eps', '-depsc');

end % scenario2


%% Scenario 3: Tuning non-linear filters
if scenario3
%% True track
% Sampling period
T = 0.1;
f = @(x) coordinatedTurnMotion(x, T); % Motion model
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(150:450) = -pi/301/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X = zeros(length(x0),K+1);
X(:,1) = x0;

% Create true track
for i=2:K+1
    % Simulate
    X(:,i) = coordinatedTurnMotion(X(:,i-1), T);
    % Set turn−rate
    X(5,i) = omega(i);
end

% Prior
mu_0 = [0; 0; 0; 0; 0];
var_0 = diag([10 10 10 5*pi/180 pi/180].^2);

% Sensor
s1 = [300; -100];
s2 = [300; -300];
h = @(x) dualBearingMeasurement(x, s1, s2); % Measurement model
R = diag([pi/180, pi/180].^2); % Measurement noise covariance

% Measurement sequence
Y = genNonLinearMeasurementSequence(X, h, R); % Generate the measurement sequence


%% a)
Q_default = diag([0, 0, 1, 0, pi/180].^2);

Q_noisy_speed = diag([0, 0, 1*100, 0, pi/180].^2);
Q_noisy_turn_rate = diag([0, 0, 1, 0, pi/180*100].^2);
Q_noisy_both = diag([0, 0, 1*100, 0, pi/180*100].^2);

Q_good_speed = diag([0, 0, 1/100, 0, pi/180].^2);
Q_good_turn_rate = diag([0, 0, 1, 0, pi/180/100].^2);
Q_good_both = diag([0, 0, 1/100, 0, pi/180/100].^2);

Q_to_test.data = {Q_default, Q_noisy_speed, Q_noisy_turn_rate, Q_noisy_both, Q_good_speed, Q_good_turn_rate, Q_good_both};
Q_to_test.names = {"Q_default", "Q_noisy_speed", "Q_noisy_turn_rate", "Q_noisy_both", "Q_good_speed", "Q_good_turn_rate", "Q_good_both"};

figure;
plotCkfEstimates(Q_to_test.data{1}, s1, s2, X, Y, mu_0, var_0, f, h, R);
xlim([-30 510]);
ylim([-460 105]);
title(Q_to_test.names{1}, 'Interpreter', 'none');
print([strcat('Images/3_a_', Q_to_test.names{1}, '.eps')], '-depsc');
legend('True state', 'CKF estimates', 'True initial position', 'Predicted initial position', 'Measurements', 'Location', 'eastoutside');

figure;
for i = 2:3
    subplot(1, 2,i-1);
    plotCkfEstimates(Q_to_test.data{i}, s1, s2, X, Y, mu_0, var_0, f, h, R);
    xlim([-30 510]);
    ylim([-460 105]);
    title(Q_to_test.names{i}, 'Interpreter', 'none');
end
print([strcat('Images/3_a_', Q_to_test.names{i-1}, '-', Q_to_test.names{i}, '.eps')], '-depsc');

figure;
plotCkfEstimates(Q_to_test.data{4}, s1, s2, X, Y, mu_0, var_0, f, h, R);
xlim([-30 510]);
ylim([-460 105]);
title(Q_to_test.names{4}, 'Interpreter', 'none');
print([strcat('Images/3_a_', Q_to_test.names{4}, '.eps')], '-depsc');


figure;
for i = 5:6
    subplot(1, 2,i-4);
    plotCkfEstimates(Q_to_test.data{i}, s1, s2, X, Y, mu_0, var_0, f, h, R);
    xlim([-30 510]);
    ylim([-460 105]);
    title(Q_to_test.names{i}, 'Interpreter', 'none');
end
print([strcat('Images/3_a_', Q_to_test.names{i-1}, '-', Q_to_test.names{i}, '.eps')], '-depsc');

figure;
plotCkfEstimates(Q_to_test.data{7}, s1, s2, X, Y, mu_0, var_0, f, h, R);
xlim([-30 510]);
ylim([-460 105]);
title(Q_to_test.names{7}, 'Interpreter', 'none');
print([strcat('Images/3_a_', Q_to_test.names{7}, '.eps')], '-depsc');

%% b) Tuned Q (X_true has constant speed )
Q_tuned = diag([0, 0, 0.1, 0, pi/180/10].^2);
[xf, Pf, xp, Pp] = nonLinearKalmanFilter(Y, mu_0, var_0, f, Q_tuned, h, R, 'CKF');

figure;
plotCkfEstimates(Q_tuned, s1, s2, X, Y, mu_0, var_0, f, h, R);
xlim([-30 510]);
ylim([-460 105]);
title("Q_tuned", 'Interpreter', 'none');
print([strcat('Images/3_b_Q_tuned.eps')], '-depsc');




%% c) Plot contours
figure;
subplot(2, 1, 1);
[xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y, mu_0, var_0, f, Q_noisy_both, h, R, 'CKF');
plot(X(1,:), X(2,:), 'r--x');
hold on;
plot(xf_ckf(1,:), xf_ckf(2,:), 'g--x');
plot(X(1,1), X(2,1), 'ro'); % Add a circle at the first point of X_true
plot(xf_ckf(1,1), xf_ckf(2,1), 'black o'); % Add a circle at the first point of X_true
N = size(Y, 2);
sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for j=1:N
    [sensor_x, sensor_y] = getPosFromMeasurement(Y(1,j), Y(2,j), s1, s2);
    % Filter out the measurements that are outside the area of interest
    if abs(sensor_x) <= 500 && abs(sensor_y) <= 500
        sensorPositions_x(j) = sensor_x;
        sensorPositions_y(j) = sensor_y;
    end
end
plot(sensorPositions_x, sensorPositions_y, '.', 'Color', [.4 .4 .4], 'MarkerSize', 8);
hold off;
xlabel('x');
ylabel('y');

xlim([-30 510]);
ylim([-460 105]);
title("Q_noisy_both", 'Interpreter', 'none');
hold on;
% Plot the 3-sigma level curve around every 5th point
for i=5:5:N
    xy_3_ckf = sigmaEllipse2D(xf_ckf(1:2,i), Pf_ckf(1:2,1:2,i), 3, 100);
    plot(xy_3_ckf(1,:), xy_3_ckf(2,:), 'g-');
end
hold off;

subplot(2, 1, 2);
%% Plot position error vs time
plot(sqrt(sum((X(1:2,2:end) - xf_ckf(1:2,:)).^2)), 'r');
xlabel('Time');
ylabel('Position error');


print('Images/3_c_noisy_both.eps', '-depsc');

figure;
subplot(2, 1, 1);
[xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y, mu_0, var_0, f, Q_good_both, h, R, 'CKF');
plot(X(1,:), X(2,:), 'r--x');
hold on;
plot(xf_ckf(1,:), xf_ckf(2,:), 'g--x');
plot(X(1,1), X(2,1), 'ro'); % Add a circle at the first point of X_true
plot(xf_ckf(1,1), xf_ckf(2,1), 'black o'); % Add a circle at the first point of X_true
N = size(Y, 2);
sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for j=1:N
    [sensor_x, sensor_y] = getPosFromMeasurement(Y(1,j), Y(2,j), s1, s2);
    % Filter out the measurements that are outside the area of interest
    if abs(sensor_x) <= 500 && abs(sensor_y) <= 500
        sensorPositions_x(j) = sensor_x;
        sensorPositions_y(j) = sensor_y;
    end
end
plot(sensorPositions_x, sensorPositions_y, '.', 'Color', [.4 .4 .4], 'MarkerSize', 8);
hold off;
xlabel('x');
ylabel('y');

xlim([-30 510]);
ylim([-460 105]);
title("Q_good_both", 'Interpreter', 'none');
hold on;
% Plot the 3-sigma level curve around every 5th point
N = size(X, 2);
for i=5:5:N
    xy_3_ckf = sigmaEllipse2D(xf_ckf(1:2,i), Pf_ckf(1:2,1:2,i), 3, 100);
    plot(xy_3_ckf(1,:), xy_3_ckf(2,:), 'g-');
end
hold off;

subplot(2, 1, 2);
%% Plot position error vs time
plot(sqrt(sum((X(1:2,2:end) - xf_ckf(1:2,:)).^2)), 'r');
xlabel('Time');
ylabel('Position error');

print('Images/3_c_good_both.eps', '-depsc');

figure;
subplot(2, 1, 1);
[xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y, mu_0, var_0, f, Q_tuned, h, R, 'CKF');
plot(X(1,:), X(2,:), 'r--x');
hold on;
plot(xf_ckf(1,:), xf_ckf(2,:), 'g--x');
plot(X(1,1), X(2,1), 'ro'); % Add a circle at the first point of X_true
plot(xf_ckf(1,1), xf_ckf(2,1), 'black o'); % Add a circle at the first point of X_true
N = size(Y, 2);
sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for j=1:N
    [sensor_x, sensor_y] = getPosFromMeasurement(Y(1,j), Y(2,j), s1, s2);
    % Filter out the measurements that are outside the area of interest
    if abs(sensor_x) <= 500 && abs(sensor_y) <= 500
        sensorPositions_x(j) = sensor_x;
        sensorPositions_y(j) = sensor_y;
    end
end
plot(sensorPositions_x, sensorPositions_y, '.', 'Color', [.4 .4 .4], 'MarkerSize', 8);
hold off;
xlabel('x');
ylabel('y');

xlim([-30 510]);
ylim([-460 105]);
title("Q_tuned", 'Interpreter', 'none');
hold on;
% Plot the 3-sigma level curve around every 5th point
N = size(X, 2);
for i=5:5:N
    xy_3_ckf = sigmaEllipse2D(xf_ckf(1:2,i), Pf_ckf(1:2,1:2,i), 3, 100);
    plot(xy_3_ckf(1,:), xy_3_ckf(2,:), 'g-');
end
hold off;

subplot(2, 1, 2);
%% Plot position error vs time
plot(sqrt(sum((X(1:2,2:end) - xf_ckf(1:2,:)).^2)), 'r');
xlabel('Time');
ylabel('Position error');

print('Images/3_c_tuned.eps', '-depsc');



end % scenario3



%% Export the source code as .txt file.
filename = fullfile('main.m');
copyfile(filename,'main.txt','f');