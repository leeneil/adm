new_size = 101;
padded_size = new_size * 5;

img0 = imread('lena30.jpg');
img_gray = rgb2gray( img0 );

img_new_size = imresize( img_gray, [new_size new_size]);

pad_size = (padded_size - new_size) / 2;

img_padded = padarray( img_new_size, [pad_size pad_size]);

mask = img_padded > 0;

save(['lena_' int2str(new_size)], 'img_padded', 'mask');




figure(1001);
imshow(img_padded);