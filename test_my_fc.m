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

delete('Images/*.eps');

% addpath('HA1');
% addpath('HA2');
% addpath('HA3');
addpath('HA4');

% rng(2); % Set random seed

scenario1 = true;
scenario2 = true;
scenario3 = true;

if scenario1

tol = 0.5;

% Number of time steps;
N = 100;

% Define prior
x_0     = [0 0 10 0 0]';
n       = length(x_0);
P_0     = diag([1 1 1 1*pi/180 1*pi/180].^2);

% Covariance
sigV = 1;
sigOmega = 1*pi/180;
G = [zeros(2,2); 1 0; 0 0; 0 1];
Q = G*diag([sigV^2 sigOmega^2])*G';

% Motion model
motionModel = @coordinatedTurnMotion;

% Random sensor position sequence
S = zeros(2,N);

% Measurement noise covariance
R = diag([10 5*pi/180].^2).*1000;

% Measurement model
measModel = @rangeBearingMeasurements;

% function handle for generating sigma points
genSigmaPoints = @sigmaPoints;


% Sample time
T = rand;

% generate state sequence
X = genNonLinearStateSequence(x_0, P_0, motionModel, T, Q, N);

% generate measurements
Y = genNonLinearMeasurementSequence(X, S, measModel, R);

% Kalman filter
[xEs, PEs, xEf, PEf, xEp, PEp] = ...
    nonLinRTSsmoother(Y, x_0, P_0, motionModel, T, Q, S, measModel, R, genSigmaPoints, 'EKF');

% plot
figure;
hold on;
plot(X(1,:),X(2,:),'b');
plot(S(:,1),S(:,1),'ko');
plot(xEs(1,:),xEs(2,:),'g');
plot(xEf(1,:),xEf(2,:),'k');
plot(xEp(1,:),xEp(2,:),'m');
hold off;
legend('True', 'sensors position' , 'Smoothed','Filtered','Predicted');
xlabel('X-position');
ylabel('Y-position');
grid on;


end % scenario1



%% Export the source code as .txt file.
filename = fullfile('main.m');
copyfile(filename,'main.txt','f');





%% TO REMOVE


function [SP,W] = sigmaPoints(x, P, type)
    % SIGMAPOINTS computes sigma points, either using unscented transform or
    % using cubature.
    %
    %Input:
    %   x           [n x 1] Prior mean
    %   P           [n x n] Prior covariance
    %
    %Output:
    %   SP          [n x 2n+1] matrix with sigma points
    %   W           [1 x 2n+1] vector with sigma point weights
    %

        switch type
            case 'UKF'

                % Dimension of state
                n = length(x);

                % Allocate memory
                SP = zeros(n,2*n+1);

                % Weights
                W = [1-n/3 repmat(1/6,[1 2*n])];

                % Matrix square root
                sqrtP = sqrtm(P);

                % Compute sigma points
                SP(:,1) = x;
                for i = 1:n
                    SP(:,i+1) = x + sqrt(1/2/W(i+1))*sqrtP(:,i);
                    SP(:,i+1+n) = x - sqrt(1/2/W(i+1+n))*sqrtP(:,i);
                end

            case 'CKF'

                % Dimension of state
                n = length(x);

                % Allocate memory
                SP = zeros(n,2*n);

                % Weights
                W = repmat(1/2/n,[1 2*n]);

                % Matrix square root
                sqrtP = sqrtm(P);

                % Compute sigma points
                for i = 1:n
                    SP(:,i) = x + sqrt(n)*sqrtP(:,i);
                    SP(:,i+n) = x - sqrt(n)*sqrtP(:,i);
                end

            otherwise
                error('Incorrect type of sigma point')
        end
    end

    function [h, H] = rangeBearingMeasurements(x, s)
    %RANGEBEARINGMEASUREMENTS calculates the range and the bearing to the
    %position given by the state vector x, from a sensor locateed in s
    %
    %Input:
    %   x           [n x 1] State vector
    %   s           [2 x 1] Sensor position
    %
    %Output:
    %   h           [2 x 1] measurement vector
    %   H           [2 x n] measurement model Jacobian
    %
    % NOTE: the measurement model assumes that in the state vector x, the first
    % two states are X-position and Y-position.

        % Range
        rng = norm(x(1:2)-s);
        % Bearing
        ber = atan2(x(2)-s(2),x(1)-s(1));
        % Measurement vector
        h = [rng;ber];

        % Measurement model Jacobian
        H = [
            (x(1)-s(1))/rng      (x(2)-s(2))/rng     0 0 0;
            -(x(2)-s(2))/(rng^2) (x(1)-s(1))/(rng^2) 0 0 0
            ];

    end

    function [f, F] = coordinatedTurnMotion(x, T)
    %COORDINATEDTURNMOTION calculates the predicted state using a coordinated
    %turn motion model, and also calculated the motion model Jacobian
    %
    %Input:
    %   x           [5 x 1] state vector
    %   T           [1 x 1] Sampling time
    %
    %Output:
    %   f           [5 x 1] predicted state
    %   F           [5 x 5] motion model Jacobian
    %
    % NOTE: the motion model assumes that the state vector x consist of the
    % following states:
    %   px          X-position
    %   py          Y-position
    %   v           velocity
    %   phi         heading
    %   omega       turn-rate

        % Velocity
        v = x(3);
        % Heading
        phi = x(4);
        % Turn-rate
        omega = x(5);

        % Predicted state
        f = x + [
            T*v*cos(phi);
            T*v*sin(phi);
            0;
            T*omega;
            0];

        % Motion model Jacobian
        F = [
            1 0 T*cos(phi) -T*v*sin(phi) 0;
            0 1 T*sin(phi) T*v*cos(phi)  0;
            0 0 1          0             0;
            0 0 0          1             T;
            0 0 0          0             1
            ];
    end

    function X = genNonLinearStateSequence(x_0, P_0, f, T, Q, N)
    %GENLINEARSTATESEQUENCE generates an N-long sequence of states using a
    %    Gaussian prior and a linear Gaussian process model
    %
    %Input:
    %   x_0         [n x 1] Prior mean
    %   P_0         [n x n] Prior covariance
    %   f           Motion model function handle
    %   T           Sampling time
    %   Q           [n x n] Process noise covariance
    %   N           [1 x 1] Number of states to generate
    %
    %Output:
    %   X           [n x N] State vector sequence
    %

        % Dimension of state vector
        n = length(x_0);

        % allocate memory
        X = zeros(n, N);

        % Generete start state
        X(:,1) = mvnrnd(x_0', P_0)';

        % Generate sequence
        for k = 2:N+1

            % generate noise vector
            q = mvnrnd(zeros(1,n), Q)';

            % Propagate through process model
            [fX, ~] = f(X(:,k-1),T);
            X(:,k) = fX + q;

        end

    end

    function Y = genNonLinearMeasurementSequence(X, S, h, R)
    %GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states
    % sequence X using a non-linear measurement model.
    %
    %Input:
    %   X           [n x N+1] State vector sequence
    %   S           [n x N] Sensor position vector sequence
    %   h           Measurement model function handle
    %   R           [m x m] Measurement noise covariance
    %
    %Output:
    %   Y           [m x N] Measurement sequence
    %

        % Parameters
        N = size(X,2);
        m = size(R,1);

        % Allocate memory
        Y = zeros(m,N-1);

        for k = 1:N-1
            % Measurement
            [hX,~] = h(X(:,k+1),S(:,k));
            % Add noise
            Y(:,k) = hX + mvnrnd(zeros(1,m), R)';

        end

    end