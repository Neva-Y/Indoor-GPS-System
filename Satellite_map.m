close all
clc

img = imread('oldeng.png');
flipdim(img,1);
img = imrotate(img, 60)
imagesc([-45 5], [-60 10],img);
hold on
[trueX, trueY] = Spherical2Azimuth(Position.latitude, Position.longitude, Position.latitude(1), Position.longitude(1), 0, 0, 6371000*pi);

plot(trueX, trueY, 'b-*','linewidth',1.5)
