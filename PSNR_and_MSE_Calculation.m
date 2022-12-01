clc
clear all
close all
read_image_1 = imread('fall_frames74.jpg');
max_element = max(read_image_1(:));
read_image_2 = imread('Noisy_Image27.jpg');
SE = ((read_image_1 - read_image_2).^2);
SE_vec = SE(:);
SE_sum = sum(SE_vec);
MSE = SE_sum./307200
%PSNR = 10 * log(255^2 ./ MSE)
PSNR = psnr(read_image_2, read_image_1)
