function plotCkfEstimates(Q_to_test, s1, s2, X, Y, mu_0, var_0, f, h, R)
    [xf_ckf, Pf_ckf, xp_ckf, Pp_ckf] = nonLinearKalmanFilter(Y, mu_0, var_0, f, Q_to_test, h, R, 'CKF');
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

end