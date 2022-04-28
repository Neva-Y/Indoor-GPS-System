function [vectorLen, x, lat, long, gyroAzimuth, azimuth, Orientationx, t] = heading_direction(file)
    load(file);
    lat = Position.latitude;
    long = Position.longitude;

    time = MagneticField.Timestamp;
    vectorLen = min([length(time), length(AngularVelocity.X), ...
        length(Orientation.X), length(MagneticField.X)]);

    Bx = MagneticField.X;
    By = MagneticField.Y;
    Bz = MagneticField.Z;
    B = [Bx By Bz];
    B = B*10^-6;

    [A,b] = magcal(B);
    B_cal = (B-b)*A;
    Bx = B_cal(:,1);
    Bx = Bx(1:vectorLen);
    By = B_cal(:,2);
    By = By(1:vectorLen);
    azimuth = (2*pi - atan2(By, Bx));
    azimuth = unwrap(azimuth) * 180/pi;
    initialHeading = azimuth(1);


    Orientationx = -Orientation.X;
    Orientationx = Orientationx(1:vectorLen);
    offset = Orientationx(1) - initialHeading;
    Orientationx = unwrap((Orientationx - offset)*pi/180);
    Orientationx = Orientationx * 180/pi;

    Wz = AngularVelocity.Z;

    % Gain for yaw due to to deflection of vertical axis due to the Earth
    % elipsoid
    yawGain = 1.1;
    Wz = Wz(1:vectorLen) * yawGain;



    yawEstimate = cumtrapz(Wz);
    gyroAzimuth = unwrap((initialHeading + yawEstimate) * pi/180);
    gyroAzimuth = gyroAzimuth * 180/pi;

    initial_min = double(minute(time(1)));
    initial_sec = double(second(time(1)));

    t = zeros(vectorLen,1);
    for i=1:length(time)
        currentmin = double(minute(time(i)));
        if currentmin == initial_min + 1
            t(i) = double(second(time(i))) + 60 - initial_sec;
        elseif currentmin == initial_min + 2
            t(i) = double(second(time(i))) + 120 - initial_sec;
        elseif currentmin == initial_min + 3
            t(i) = double(second(time(i))) + 180 - initial_sec;
        elseif currentmin == initial_min + 4
            t(i) = double(second(time(i))) + 240 - initial_sec;    
        else
            t(i) = double(second(time(i)))-initial_sec;
        end
    end

    %Sampling rate
    Fs = 50;
    dt = 1/Fs;

    % Q is the covariance of the gaussian noise for the modelling errors
    % R is the covariance of the gaussian noise for the measurement noise
    Q = 0.0003;
    R = 1.0585;

    %Initialise Covariance
    P = 0;

    %Initial state and state space model
    z = azimuth;
    x = zeros(vectorLen,1);   
    u = Wz;
    x(1) = initialHeading;
    F = 1;
    G = dt;
    H = 1;

    t = t(1:vectorLen);

    for i = 2:vectorLen
        % Prediction for state vector and covariance
        x(i) = F * x(i-1) + G * u(i-1);
        P = F * P * F.' + Q;
        % Computing the Kalman gain, closer to 1 means a greater trust in
        % measurements
        K = P * H.' * inv(H * P * H.' + R);

        % Correction based on observation
        x(i) = x(i) + K * (z(i) - H*x(i));
        P = P - K * H * P;
    end


    % fig1 = figure();
    % plot(t,BBx, 'LineWidth', 1)
    % hold on
    % grid on
    % plot(t,BBy, 'LineWidth', 1)
    % plot(t,BBz, 'LineWidth', 1)
    % xlabel('Time (sec)');
    % ylabel('Magnetic Field (T)');
    % legend('x','y','z');
    % title("Magnetic Field measured");
    % 
    % fig2 = figure();
    % plot(t,W(:,1), 'LineWidth', 1)
    % hold on
    % grid on
    % plot(t,W(:,2), 'LineWidth', 1)
    % plot(t,W(:,3), 'LineWidth', 1)
    % xlabel('Time (sec)');
    % ylabel('Angular Velocity (deg/s)');
    % legend('x','y','z');
    % title("Angular Velocity measured");
    % 
    % 
    % fig3 = figure();
    % plot(t, azimuth, 'LineWidth', 1);
    % xlabel('Time (sec)');
    % ylabel('Heading Angle (deg)');
    % title("Heading Direction Relative to North using Magnetometer");
    % 
    % fig4 = figure();
    % plot(t, gyroAzimuth, 'LineWidth', 1);
    % xlabel('Time (sec)');
    % ylabel('Heading Angle (deg)');
    % title("Heading Direction Relative to North using Gyro");
    % 
    % fig5 = figure();
    % plot(t, Orientationy, 'LineWidth', 1);
    % xlabel('Time (sec)');
    % ylabel('Heading Angle (deg)');
    % title("True Heading Direction");

%     fig6 = figure();
%     plot(t, gyroAzimuth, 'LineWidth', 1);
%     hold on 
%     plot(t, azimuth, 'LineWidth', 1);
%     plot(t, Orientationx, 'LineWidth', 1);
%     plot(t, x, 'k','linewidth',2);
%     title("Heading Direction");
%     xlabel('Time (sec)');
%     ylabel('Heading Angle (deg)');
%     legend('Gyro','Magnetometer','True Heading', 'Kalman Prediction')
end