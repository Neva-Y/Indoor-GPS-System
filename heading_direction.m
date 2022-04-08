clc
clear all
close all
load('/Users/yang/MATLAB-Drive/MobileSensorData/gyro_mag_spin.mat')

time = MagneticField.Timestamp
Bx = MagneticField.X;
By = MagneticField.Y;
Bz = MagneticField.Z;
B = [Bx By Bz];

[A,b,expmfs] = magcal(B);
B_cal = (B-b)*A;

Wx = AngularVelocity.X;
Wy = AngularVelocity.Y;
Wz = AngularVelocity.Z;
W = [Wx Wy Wz];


initial_min = double(minute(time(1)));
initial_sec = double(second(time(1)));

for i=1:length(time)
    if double(minute(time(i)))>initial_min
        t(i) = double(second(time(i))) + 60 - initial_sec;
    else
        t(i) = double(second(time(i)))-initial_sec;
    end
end

BBx = B_cal(:,1);
BBy = B_cal(:,2);
BBz = B_cal(:,3);

fig1 = figure();
plot(t,B_cal(:,1), 'LineWidth', 1)
hold on
grid on
plot(t,B_cal(:,2), 'LineWidth', 1)
plot(t,B_cal(:,3), 'LineWidth', 1)
xlabel('Time (sec)');
ylabel('Magnetic Field (uT)');
legend('x','y','z');
title("Magnetic Field measured");

fig2 = figure();
plot(t,W([1:870],1), 'LineWidth', 1)
hold on
grid on
plot(t,W([1:870],2), 'LineWidth', 1)
plot(t,W([1:870],3), 'LineWidth', 1)
xlabel('Time (sec)');
ylabel('Angular Velocity (rad/s)');
legend('x','y','z');
title("Angular Velocity measured");


azimuth = atan2(BBx, BBy) * 180.0/pi;
fig3 = figure();
plot(t, azimuth, 'LineWidth', 1);
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');
title("Heading Direction Relative to North");


