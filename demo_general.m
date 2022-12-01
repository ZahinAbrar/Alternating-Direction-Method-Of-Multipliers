clear all
close all
clc

% read test image
% this has been commented to reduce the size of the image
%  z = im2double(imread('frames100.jpg'));      %change
% imageArray_resize = imageArray([601:1000],[601:1000])

read_image = imread('fall_frames149.jpg')
% imageArray_resize = read_image([601:1000],[601:1000]); % when I read the whole image the imshow do not display full picture shows 67% of it
z = im2double(read_image)


%construct A matrix, deblurring as an example
dim = size(z);

h = fspecial('gaussian',7,.0000001);       % I am tuning this one, originally it was h = fspecial('gaussian',[9,9],1); Quick Result 7,.0000001, [480,640],1
A = @(z,trans_flag) afun(z,trans_flag,h,dim);

%reset random number generator
rng(0);

%set noies level
noise_level = 0;

%calculate observed image
y = A(z(:),'transp') + noise_level*randn(prod(dim),1);   % observing in vector fromat
y = proj(y,[0,1]);
y = reshape(y,dim);                                      % reshaping the image format


% 
% %initialize variables
N           = dim(1)*dim(2);    %number of pixels in the image
v           = 0.5*ones(dim);
x           = v;
u           = zeros(dim);
residual    = inf;
% 
% % In first iteration we need to uncomment these three lines below
%%%
% x_old = x;  % or x_update for the first iteration
% v_old = v;  % or v_update for the second iteration
u_old = u;  % or u_update for the third iteration
%%% 
 
% %optional parameters

opts.rho     = 1;
opts.gamma   = 1;
opts.max_itr = 20;
opts.print   = true;
 
% % % % % % the four lines below will be uncommented from the second iteration
% % % % % % and go on till the last iteration


% load('v_update')        
% load('u_update')
% v = v_update;
% u = u_update;
%  
% %%% main routine
% 
% x_update = Inversion_Step(y,A,opts,v,u);
% Noisy_Output = x_update + u;
% imwrite(Noisy_Output ,'Noisy_Image.jpg')
% save('x_update.mat', 'x_update')
% save('Noisy_Output.mat', 'Noisy_Output') % because this noisy op will go to the denoiser

%%%

%%%Now we need to read the denoised image for the third step of the ADMM %%%

%% Second step v_update
% 
load('Noisy_Output')
v_update = denoise(Noisy_Output);
save('v_update.mat','v_update')

%% End Second step
 
%% Third step u_update
% 
% load('v_update');
% load('x_update');
% % load('Noisy_Output');
% x_update = x_update;
% v_update = v_update;
% % % % this and below line will be commented in the first iteration
% load('u_update') 
% u_old = u_update;
% % % %
% u_update = final(u_old, x_update, v_update);
% save('u_update.mat','u_update')


%% End Third step








 
%% A gentle Guideline for Running the code %%
% the one problem I faced I couldn't manage to run the code for denoiser
% automatically in python, that's why this code need to be run manually
% with a great care. There are three section of the code inversing
% (x_updating), denoising(v_updating) and u_updating. First we run the
% inversing, then we save the noisy image in a folder and manually pass
% through our autoneocoder denoiser. Then we save the denoised output and
% read it to perform the v-updating while commenting out the inversion
% block. In both of the process we save the updated x and updated v
% manually. After getting the updated_v we commneted the denoising a.k.a
% v_updating block. Finally we go for the u_updating. That's how the
% process go on. 
% at the very first iteration v_old and u_old are just some predefined
% value, so be careful while commenting out different portion. But once an
% iteration is done saved "updated_x, updated_v and updated_u" should be
% used to get the algorithm-wise result.


% %display
%later Needed
% PSNR_output = psnr(Noisy_output,z);
% fprintf('\nPSNR = %3.2f dB \n', PSNR_output);
% 
% figure;
% subplot(121);
% imshow(y);
% title('Input');
% 
% subplot(122);
% imshow(Noisy_output);
% tt = sprintf('PSNR = %3.2f dB', PSNR_output);
% title(tt);
