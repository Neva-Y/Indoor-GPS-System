close all
clear all
clc
load("Data/oldeng2.mat");
    lat = Position.latitude;
    long = Position.longitude;

    time = MagneticField.Timestamp;
    vectorLen = min([length(time), length(AngularVelocity.X), ...
        length(Orientation.X), length(MagneticField.X)]);

    Bx = MagneticField.X*10^-6;
    By = MagneticField.Y*10^-6;
    Bz = MagneticField.Z*10^-6;
    B = [Bx By Bz];

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
    SSE_low = 100000;
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
    yawGain=1.05
    Wz = Wz(1:vectorLen) * yawGain;
        yawEstimate = cumtrapz(Wz);
        gyroAzimuth = unwrap((initialHeading + yawEstimate) * pi/180);
        gyroAzimuth = gyroAzimuth * 180/pi;
    fig6 = figure();
    
    plot(t, gyroAzimuth, 'LineWidth', 1);
    hold on
    plot(t, Orientationx, 'LineWidth', 1);
    title("Heading Direction");
    xlabel('Time (sec)');
    ylabel('Heading Angle (deg)');
    legend('Gyro with Gain k = 1.05','True Heading')


  