clc
clear all
close all
load('Data/heading_north.mat');

time = MagneticField.Timestamp;
vectorLen = length(time);

Orientationx = -Orientation.X;
Orientationx = Orientationx(1:vectorLen);
Orientationy = Orientation.Y;
Orientationy = Orientationy(1:vectorLen);
Orientationz = Orientation.Z;
Orientationz = Orientationz(1:vectorLen);


Bx = MagneticField.X;
By = MagneticField.Y;
Bz = MagneticField.Z;
B = [Bx By Bz];
B = B * 10^-6;

[A,b,expmfs] = magcal(B);
B_cal = (B-b)*A;
BBx = B_cal(:,1);
BBx = BBx(1:vectorLen);
BBy = B_cal(:,2);
BBy = BBy(1:vectorLen);
BBz = B_cal(:,3);
BBz = BBz(1:vectorLen);
azimuth = atan2(BBx, BBy) * 180.0/pi;
initialHeading = azimuth(1);

Wx = AngularVelocity.X;
Wx = Wx(1:vectorLen);
Wy = AngularVelocity.Y;
Wy = Wy(1:vectorLen);
Wz = AngularVelocity.Z;
Wz = Wz(1:vectorLen);
W = [Wx Wy Wz];
W = W * 180.0/pi;
yawEstimate = cumtrapz(Wz);
gyroAzimuth = initialHeading + yawEstimate;


initial_min = double(minute(time(1)));
initial_sec = double(second(time(1)));

for i=1:length(time)
    if double(minute(time(i)))>initial_min
        t(i) = double(second(time(i))) + 60 - initial_sec;
    else
        t(i) = double(second(time(i)))-initial_sec;
    end
end



fig1 = figure();
plot(t,B_cal(:,1), 'LineWidth', 1)
hold on
grid on
plot(t,B_cal(:,2), 'LineWidth', 1)
plot(t,B_cal(:,3), 'LineWidth', 1)
xlabel('Time (sec)');
ylabel('Magnetic Field (T)');
legend('x','y','z');
title("Magnetic Field measured");

fig2 = figure();
plot(t,W(:,1), 'LineWidth', 1)
hold on
grid on
plot(t,W(:,2), 'LineWidth', 1)
plot(t,W(:,3), 'LineWidth', 1)
xlabel('Time (sec)');
ylabel('Angular Velocity (deg/s)');
legend('x','y','z');
title("Angular Velocity measured");


fig3 = figure();
plot(t, azimuth, 'LineWidth', 1);
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');
title("Heading Direction Relative to North using Magnetometer");

fig4 = figure();
plot(t, gyroAzimuth, 'LineWidth', 1);
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');
title("Heading Direction Relative to North using Gyro");

fig5 = figure();
plot(t, Orientationx, 'LineWidth', 1);
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');
title("True Heading Direction");

fig5= figure();
plot(t, gyroAzimuth, 'LineWidth', 1);
hold on 
plot(t, azimuth, 'LineWidth', 1);
plot(t, Orientationx, 'LineWidth', 1);
title("Heading Direction");
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');
legend('Gyro','Magnetometer','True')