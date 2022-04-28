close all
clc

img = imread('oldeng2.png');
img = imrotate(img,0,'bilinear','crop');
%imagesc([-37.80005 -37.7988], [144.9609 144.9626],img);
imagesc([-77 70], [-39 98], flipud(img));


hold on
[trueX, trueY] = Spherical2Azimuth(Position.latitude, Position.longitude, Position.latitude(1), Position.longitude(1), 0, 0, 6371000*pi);

plot(trueX, trueY, 'b-*','linewidth',1.5);
set(gca,'ydir','normal');
