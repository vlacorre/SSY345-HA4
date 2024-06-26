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

clc; clear all; close all;

delete('Images/*.eps');

addpath('HA1');
addpath('HA2');
addpath('HA3');
addpath('HA4');

rng(1); % Set random seed

scenario1 = false;
scenario2 = false;
scenario3 = true;

if scenario1
%% a)

%% True track
% Sampling period
T = 0.1;
% Length of time sequence
K = 600;
% Allocate memory
omega = zeros(1,K+1);
% Turn rate
omega(150:450) = -pi/301/T;
% Initial state
x0 = [0 0 20 0 omega(1)]';
% Allocate memory
X_true = zeros(length(x0),K+1);
X_true(:,1) = x0;

% Create true track
for i=2:K+1
    % Simulate
    X_true(:,i) = coordinatedTurnMotion(X_true(:,i-1), T);
    % Set turn−rate
    X_true(5,i) = omega(i);
end
% Process model
f = @(x, T) coordinatedTurnMotion(x, T); % Motion model
Q_tuned = diag([0 0 0.1 0 pi/1800].^2);

% Prior
x_0 = zeros(5,1);
P_0 = diag([10 10 10 5*pi/180 pi/180].^2);

% Measurement model
s1 = [300; -100];
s2 = [300; -300];
S = [s1(1) s2(1) 0 0 0;s1(2) s2(2) 0 0 0];
h = @(x) dualBearingMeasurement(x, s1, s2); % Measurement model
sigma_phi_1 = pi/180;
sigma_phi_2 = pi/180;
R = diag([sigma_phi_1 sigma_phi_2].^2); % Measurement noise covariance

Y = genNonLinearMeasurementSequence(X_true, h, R);
h = @(x, S) dualBearingMeasurement(x, S(:,1), S(:,2)); % Measurement model

% Positions corresponding to the measurements
N = size(Y, 2);
sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for j=1:N
    [sensor_x, sensor_y] = getPosFromMeasurement(Y(1,j), Y(2,j), s1, s2);
    sensorPositions_x(j) = sensor_x;
    sensorPositions_y(j) = sensor_y;
end
sigmaPoints = @sigmaPoints;
% Filtered & smoothed position trajectories
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, T, Q_tuned, S, h, R, sigmaPoints, 'CKF');


% Plot
figure;
grid on;
hold on;
plot(sensorPositions_x, sensorPositions_y, '.', 'Color', [.4 .4 .4], 'MarkerSize', 8);
plot(xf(1,:), xf(2,:), 'r');
plot(xs(1,:), xs(2,:), 'b');
% Plot the 3-sigma level curve around every 5th point
for i=5:5:N
    xy_3_f = sigmaEllipse2D(xf(1:2,i), Pf(1:2,1:2,i), 3, 100);
    xy_3_s = sigmaEllipse2D(xs(1:2,i), Ps(1:2,1:2,i), 3, 100);
    plot(xy_3_f(1,:), xy_3_f(2,:), 'k-', 'LineWidth', 1);
    plot(xy_3_s(1,:), xy_3_s(2,:), 'c-', 'LineWidth', 1);
end
hold off;
xlim([-10 500]);
ylim([-450 100]);
xlabel('x');
ylabel('y');
legend('Measurements', 'Filtered', 'Smoothed', 'Filtered covariance (3\sigma)', 'Smoothed covariance (3\sigma)')

print('Images/1_a_smoothed.eps', '-depsc');

%% b) Outlier
% Modify measurement at k = 300
Y(:,300) = Y(:,300) + mvnrnd([0 0], R*1000)';

% Positions corresponding to the measurements
N = size(Y, 2);
sensorPositions_x = zeros(N, 1);
sensorPositions_y = zeros(N, 1);
for j=1:N
    [sensor_x, sensor_y] = getPosFromMeasurement(Y(1,j), Y(2,j), s1, s2);
    sensorPositions_x(j) = sensor_x;
    sensorPositions_y(j) = sensor_y;
end
sigmaPoints = @sigmaPoints;
% Filtered & smoothed position trajectories
[xs, Ps, xf, Pf, xp, Pp] = nonLinRTSsmoother(Y, x_0, P_0, f, T, Q_tuned, S, h, R, sigmaPoints, 'CKF');


% Plot
figure;
grid on;
hold on;
plot(sensorPositions_x, sensorPositions_y, '.', 'Color', [.4 .4 .4], 'MarkerSize', 8);
plot(xf(1,:), xf(2,:), 'r');
plot(xs(1,:), xs(2,:), 'b');
% Plot the 3-sigma level curve around every 5th point
for i=5:5:N
    xy_3_f = sigmaEllipse2D(xf(1:2,i), Pf(1:2,1:2,i), 3, 100);
    xy_3_s = sigmaEllipse2D(xs(1:2,i), Ps(1:2,1:2,i), 3, 100);
    plot(xy_3_f(1,:), xy_3_f(2,:), 'k-', 'LineWidth', 1);
    plot(xy_3_s(1,:), xy_3_s(2,:), 'c-', 'LineWidth', 1);
end
% Add red circle around modified measurement
plot(sensorPositions_x(300), sensorPositions_y(300), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
xlim([-10 500]);
ylim([-450 100]);
xlabel('x');
ylabel('y');
legend('Measurements', 'Filtered', 'Smoothed', 'Filtered covariance (3\sigma)', 'Smoothed covariance (3\sigma)')

print('Images/1_b_outlier.eps', '-depsc');

end % scenario1

if scenario2
%% a) Resampling

% Models
A = 1;
f = @(x, A) (A*x);
Q = 1.5;

% Measurement model
H = 1;
h = @(x, H) (H*x);
R = 3;
m = size(H,1);

% Prior
x_0 = 2;
P_0 = 8;
n = size(x_0,1);

% Number of time steps
K = 30;

% Generate state and measurement sequences
X = zeros(n,K);
Y = zeros(m,K);

q = mvnrnd(zeros(1,n), Q, K)';
r = mvnrnd(zeros(1,m), R, K)';
xk = x_0;
for k = 1:K
    xk = f(xk,A) + q(:,k);
    X(:,k) = xk;

    Y(:,k) = h(xk, H) + r(:,k);
end

% Run Kalman filter
[Xf, Pf] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

% PF with poor choice of N
N = 100;
proc_f = @(X_kmin1) (f(X_kmin1, A));
meas_h = @(X_k) (h(X_k, H));
plotFunc = @(~) (0); % Define dummy function that does nothing

[xfp, Pfp] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, false, plotFunc);

[xfpr, Pfpr] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, true, plotFunc);

mse_kf = mean((Xf(1,:) - X(1,:)).^2);
mse_pf = mean((xfp(1,:) - X(1,:)).^2);
mse_pf_resampling = mean((xfpr(1,:) - X(1,:)).^2);
fprintf('N = %d Mean Square Error (KF): %.4f\n', N, mse_kf);
fprintf('N = %d Mean Square Error (PF without resampling): %.4f\n', N, mse_pf);
fprintf('N = %d Mean Square Error (PF with resampling): %.4f\n', N, mse_pf_resampling);


% PF with good choice of N
N = 1000;
proc_f = @(X_kmin1) (f(X_kmin1, A));
meas_h = @(X_k) (h(X_k, H));
plotFunc = @plotPostPdf;
[xfp, Pfp] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, false, plotFunc);

plotFunc = @(~) (0); % Define dummy function that does nothing
[xfpr, Pfpr] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, true, plotFunc);

figure;
grid on;
plot(X(1,:), 'r');
hold on;
errorbar(Xf(1,:),Pf(1,:), 'k', 'LineWidth', 1.5);
errorbar(xfp(1,:),Pfp(1,:), 'b', 'LineWidth', 1.5);
errorbar(xfpr(1,:),Pfpr(1,:), 'm--', 'LineWidth', 1.5);
plot(Y(1,:), 'x', 'Color', [.4 .4 .4], 'MarkerSize', 10);
hold off;
xlabel('k');
ylabel('x');
legend('True state', 'KF', 'PF', 'PF with resampling', 'Measurements');

print('Images/2_a_resampling.eps', '-depsc');

mse_kf = mean((Xf(1,:) - X(1,:)).^2);
mse_pf = mean((xfp(1,:) - X(1,:)).^2);
mse_pf_resampling = mean((xfpr(1,:) - X(1,:)).^2);

fprintf('N = %d Mean Square Error (KF): %.4f\n', N, mse_kf);
fprintf('N = %d Mean Square Error (PF without resampling): %.4f\n', N, mse_pf);
fprintf('N = %d Mean Square Error (PF with resampling): %.4f\n', N, mse_pf_resampling);



%% b) Incorrect prior

% Models
A = 1;
f = @(x, A) (A*x);
Q = 1.5;

% Measurement model
H = 1;
h = @(x, H) (H*x);
R = 3;
m = size(H,1);

% Prior
x_0 = 2;
P_0 = 8;
n = size(x_0,1);

% Number of time steps
K = 30;

% Generate state and measurement sequences
X = zeros(n,K);
Y = zeros(m,K);

q = mvnrnd(zeros(1,n), Q, K)';
r = mvnrnd(zeros(1,m), R, K)';
xk = x_0;
for k = 1:K
    xk = f(xk,A) + q(:,k);
    X(:,k) = xk;

    Y(:,k) = h(xk, H) + r(:,k);
end

% Incorrect prior
x_0 = -20;
P_0 = 2;

% Run Kalman filter
[Xf, Pf] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

N = 1000;
meas_h = @(X_k) (h(X_k, H));
plotFunc = @(~) (0); % Define dummy function that does nothing

[xfp, Pfp] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, false, plotFunc);

[xfpr, Pfpr] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, true, plotFunc);


figure;
grid on;
plot(X(1,:), 'r');
hold on;
errorbar(Xf(1,:),Pf(1,:), 'k', 'LineWidth', 1.5);
errorbar(xfp(1,:),Pfp(1,:), 'b', 'LineWidth', 1.5);
errorbar(xfpr(1,:),Pfpr(1,:), 'm--', 'LineWidth', 1.5);
plot(Y(1,:), 'x', 'Color', [.4 .4 .4], 'MarkerSize', 10);
hold off;
xlabel('k');
ylabel('x');
legend('True state', 'KF', 'PF', 'PF with resampling', 'Measurements');

print('Images/2_b_incorrect_prior.eps', '-depsc');

% Correct prior
x_0 = 2;
P_0 = 8;

%% c) Illustrate the particle trajectories for a PF without resampling

N = 100;
plotFunc = @plotPartTrajs;
figure;
[xfp, Pfp] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, false, plotFunc);

grid on;
plot(X(1,:), 'r');
hold on;
errorbar(Xf(1,:),Pf(1,:), 'k', 'LineWidth', 1.5);
errorbar(xfp(1,:),Pfp(1,:), 'b', 'LineWidth', 1.5);
errorbar(xfpr(1,:),Pfpr(1,:), 'm--', 'LineWidth', 1.5);
plot(Y(1,:), 'x', 'Color', [.4 .4 .4], 'MarkerSize', 10);
hold off;
xlabel('k');
ylabel('x');
legend('True state', 'KF', 'PF', 'PF with resampling', 'Measurements');

print('Images/2_c_plotPartTrajs.eps', '-depsc');

%% d) Illustrate the particle trajectories for a PF with resampling

N = 100;
plotFunc = @plotPartTrajs;
figure;
[xfp, Pfp] = pfFilter(x_0, P_0, Y, proc_f, Q, meas_h, R, N, true, plotFunc);

grid on;
plot(X(1,:), 'r');
hold on;
errorbar(Xf(1,:),Pf(1,:), 'k', 'LineWidth', 1.5);
errorbar(xfp(1,:),Pfp(1,:), 'b', 'LineWidth', 1.5);
errorbar(xfpr(1,:),Pfpr(1,:), 'm--', 'LineWidth', 1.5);
plot(Y(1,:), 'x', 'Color', [.4 .4 .4], 'MarkerSize', 10);
hold off;
xlabel('k');
ylabel('x');
legend('True state', 'KF', 'PF', 'PF with resampling', 'Measurements');

print('Images/2_d_plotPartTrajs_resampling.eps', '-depsc');

end % scenario2

if scenario3
%% a) True state
% MapProblemGetPoint
%% b)
load Xk;

% Compute velocity from positions
V = zeros(2, size(Xk,2));
for i=2:size(Xk,2)
    V(:,i) = Xk(1:2,i) - Xk(1:2,i-1);
end

% Velocity measurements
sigma_v = 0.001;
Y = V + mvnrnd([0 0], sigma_v*eye(2), size(V,2))';
figure;
plot(Y(1,:), Y(2,:), '.', 'Color', [.4 .4 .4], 'MarkerSize', 40);
xlabel('x');
ylabel('y');
print('Images/3_b_velocity_measurements.eps', '-depsc');

%% d)
% Constant velocity motion model
T = 1; % Empirically chosen
A = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
f = @(x) A*x;
sigma_pos = 0.1;
sigma_speed = 1;
Q = diag([sigma_pos sigma_pos sigma_speed sigma_speed].^2);

% Measurement model
H = [0 0 1 0; 0 0 0 1];
meas_h = @(x) H*x;
sigma_meas_speed = 0.1;
R = diag([sigma_meas_speed sigma_meas_speed].^2);

% Prior
x_0 = [Xk(:,1); 0; 0];
P_0 = diag([0 0 0.1 0.1].^2);

N = 1000;
plotFunc = @plotPartTrajs;
plotFunc = @plotPostPdf;
plotFunc = @(~) (0); % Define dummy function that does nothing
[xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, f, Q, meas_h, R, N, true, plotFunc);
% [xfp_wd, Pfp_wd, Xp_wd, Wp_wd] = pfFilter_wall_detection(x_0, P_0, Y, f, Q, meas_h, R, N, true, plotFunc);
% [xfp, Pfp] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

figure;
grid on;
plot(Xk(1,:), Xk(2,:), 'rx-');
hold on;
plot(xfp(1,:), xfp(2,:), 'bx-');
% plot(xfp_wd(1,:), xfp_wd(2,:), 'gx-');

% plot particles around some posteriors
part_group = [];
K = size(Y, 2);
for k=1:2:K
    % Sort Xp by weight
    [~, I] = sort(Wp(:,k), 'descend');
    for i=1:3
        part_group = [part_group; Xp(1,I(i),k) Xp(2,I(i),k)];
    end
end
plot(part_group(:,1), part_group(:,2), 'm*', 'MarkerSize', 12);

% plot particles around some posteriors
% part_group = [];
% K = size(Y, 2);
% for k=1:2:K
%     % Sort Xp by weight
%     [~, I] = sort(Wp_wd(:,k), 'descend');
%     for i=1:3
%         part_group = [part_group; Xp_wd(1,I(i),k) Xp_wd(2,I(i),k)];
%     end
% end
% plot(part_group(:,1), part_group(:,2), 'k.', 'MarkerSize', 30);

plot([1+1i 1+9*1i 5+9*1i])
plot([7+9*1i 11+9*1i 11+1i 7+1i]);plot([5+1i 1+1i])
plot([2+5.2*1i 2+8.3*1i 4+8.3*1i 4+5.2*1i 2+5.2*1i])%House 1
plot([2+3.7*1i 2+4.4*1i 4+4.4*1i 4+3.7*1i 2+3.7*1i])%House 2
plot([2+2*1i 2+3.2*1i 4+3.2*1i 4+2*1i 2+2*1i])%House 3
plot([5+1i 5+2.2*1i 7+2.2*1i 7+1i])%House 4
plot([5+2.8*1i 5+5.5*1i 7+5.5*1i 7+2.8*1i 5+2.8*1i])%House 5
plot([5+6.2*1i 5+9*1i]);plot([7+9*1i 7+6.2*1i 5+6.2*1i])%House 6
plot([8+4.6*1i 8+8.4*1i 10+8.4*1i 10+4.6*1i 8+4.6*1i])%House 7
plot([8+2.4*1i 8+4*1i 10+4*1i 10+2.4*1i 8+2.4*1i])%House 8
plot([8+1.7*1i 8+1.8*1i 10+1.8*1i 10+1.7*1i 8+1.7*1i])%House 9

axis([0.8 11.2 0.8 9.2])

hold off;
xlabel('x');
ylabel('y');
legend('True state', 'PF', 'PF particles');
print('Images/3_d_PF_perfect_prior.eps', '-depsc');

%% e) No prior


% Prior
x_0 = [6; 5; 0; 0];
P_0 = diag([5 5 0.1 0.1].^2);

N = 10000;
plotFunc = @plotPartTrajs;
plotFunc = @(~) (0); % Define dummy function that does nothing
[xfp, Pfp, Xp, Wp] = pfFilter(x_0, P_0, Y, f, Q, meas_h, R, N, true, plotFunc);
[xfp_wd, Pfp_wd, Xp_wd, Wp_wd] = pfFilter_wall_detection_customP(x_0, P_0, Y, f, Q, meas_h, R, N, true, plotFunc);
% [xfp, Pfp] = kalmanFilter(Y, x_0, P_0, A, Q, H, R);

figure;
grid on;
plot(Xk(1,:), Xk(2,:), 'rx-');
hold on;
plot(xfp(1,:), xfp(2,:), 'bx-');
plot(xfp_wd(1,:), xfp_wd(2,:), 'gx-');

plot([1+1i 1+9*1i 5+9*1i])
plot([7+9*1i 11+9*1i 11+1i 7+1i]);plot([5+1i 1+1i])
plot([2+5.2*1i 2+8.3*1i 4+8.3*1i 4+5.2*1i 2+5.2*1i])%House 1
plot([2+3.7*1i 2+4.4*1i 4+4.4*1i 4+3.7*1i 2+3.7*1i])%House 2
plot([2+2*1i 2+3.2*1i 4+3.2*1i 4+2*1i 2+2*1i])%House 3
plot([5+1i 5+2.2*1i 7+2.2*1i 7+1i])%House 4
plot([5+2.8*1i 5+5.5*1i 7+5.5*1i 7+2.8*1i 5+2.8*1i])%House 5
plot([5+6.2*1i 5+9*1i]);plot([7+9*1i 7+6.2*1i 5+6.2*1i])%House 6
plot([8+4.6*1i 8+8.4*1i 10+8.4*1i 10+4.6*1i 8+4.6*1i])%House 7
plot([8+2.4*1i 8+4*1i 10+4*1i 10+2.4*1i 8+2.4*1i])%House 8
plot([8+1.7*1i 8+1.8*1i 10+1.8*1i 10+1.7*1i 8+1.7*1i])%House 9

axis([0.8 11.2 0.8 9.2])

% plot particles around each posterior
K = size(Y, 2);
for k=1:5:K
    for i=1:50
        plot(Xp_wd(1,i,k), Xp_wd(2,i,k), '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 30);
    end
end
hold off;
xlabel('x');
ylabel('y');
legend('True state', 'PF', 'PF with wall detection', 'PF particles', 'PF particles with wall detection');
print('Images/3_d_PF_no_prior.eps', '-depsc');


end % scenario3

%% Export the source code as .txt file.
filename = fullfile('main.m');
copyfile(filename,'main.txt','f');

