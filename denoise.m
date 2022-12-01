function v_update = denoise(Noisy_Output)
% Denoise has been run manually
read_denoised_image = imread('Denoised.bmp')
denoi_read = im2double(read_denoised_image)
v_update = denoi_read
end
