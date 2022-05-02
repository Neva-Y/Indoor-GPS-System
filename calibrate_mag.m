load('Data/calibrate-mag')

Bx = MagneticField.X*10^-6;
By = MagneticField.Y*10^-6;
Bz = MagneticField.Z*10^-6;
B = [Bx By Bz];

[A,b] = magcal(B);
B_cal = (B-b)*A;

figure()
plot3(Bx(:),By(:),Bz(:),'LineStyle','none','Marker','X','MarkerSize',8)
hold on
grid(gca,'on')
plot3(B_cal(:,1),B_cal(:,2),B_cal(:,3),'LineStyle','none','Marker', ...
            'o','MarkerSize',8,'MarkerFaceColor','r') 
axis equal
xlabel('T')
ylabel('T')
zlabel('T')
legend('Uncalibrated Samples', 'Calibrated Samples','Location', 'southoutside')
title("Uncalibrated vs Calibrated" + newline + "Magnetometer Measurements")
hold off