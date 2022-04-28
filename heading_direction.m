load('Data/oldeng2.mat');


time = MagneticField.Timestamp;
vectorLen = min([length(time), length(AngularVelocity.X), ...
    length(Orientation.X), length(MagneticField.X)]);

Bx = MagneticField.X;
By = MagneticField.Y;
Bz = MagneticField.Z;
B = [Bx By Bz];
B = B*10^-6;

[A,b,expmfs] = magcal(B);
B_cal = (B-b)*A;
BBx = B_cal(:,1);
BBx = BBx(1:vectorLen);
BBy = B_cal(:,2);
BBy = BBy(1:vectorLen);
BBz = B_cal(:,3);
BBz = BBz(1:vectorLen);
azimuth = (2*pi - atan2(BBy, BBx));
azimuth = unwrap(azimuth) * 180/pi;
initialHeading = azimuth(1);


Orientationx = -Orientation.X;
Orientationx = Orientationx(1:vectorLen);
Orientationy = Orientation.Y;
Orientationy = Orientationy(1:vectorLen);
Orientationz = Orientation.Z;
Orientationz = Orientationz(1:vectorLen);
offset = Orientationx(1) - initialHeading;
Orientationx = Orientationx - offset;


Wx = AngularVelocity.X;
Wx = Wx(1:vectorLen);
Wy = AngularVelocity.Y;
Wy = Wy(1:vectorLen);
Wz = AngularVelocity.Z;
Wz = Wz(1:vectorLen);
W = [Wx Wy Wz];
W = W * 180.0/pi;

% Gain for yaw due to to deflection of vertical axis due to the Earth
% elipsoid
yawGain = 1;

yawEstimate = cumtrapz(Wz) * yawGain;
gyroAzimuth = initialHeading + yawEstimate;


initial_min = double(minute(time(1)));
initial_sec = double(second(time(1)));

t = zeros(vectorLen);
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

% Noise and Measurement Error
Q = 1;
R = 4000;

%Initialise Covariance
P = 1;

%Initial state
z = azimuth;
x = zeros(vectorLen);   
u = Wz;
x(1) = initialHeading;
F = 1;
G = dt;
H = 1;


t = t(1:vectorLen);

for i = 2:vectorLen
    x(i) = F * x(i-1) + G * u(i-1);
    P = F * P * F' + Q;
    K = P * H' * inv(H * P * H' + R);
    
    
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

fig6 = figure();
plot(t, gyroAzimuth, 'LineWidth', 1);
hold on 
plot(t, azimuth, 'LineWidth', 1);
plot(t, unwrap(Orientationx), 'LineWidth', 1);
plot(t, x, 'LineWidth', 3);
title("Heading Direction");
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');
legend('Gyro','Magnetometer','True Heading', 'Kalman Prediction')