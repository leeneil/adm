load('lena_101.mat');

img = im2double( img_padded );

fimg = abs( fft2( img ) );




figure(1001);
imshow(img);

figure(1002);
imagesc( log10( fftshift(fimg.^2) ) );
axis image;
colormap(jet);
colorbar;

rng(5566);

x_init = support_constraint( rand(size(img)), mask );
y_init = yconstraint( rand(size(img)), fimg );

figure(1011);
imshow(x_init);

figure(1012);
imshow(y_init);

n_iter = 500;

x = x_init;
y = y_init;
lambda = zeros(size(x));
pi0 = zeros(size(x));

beta0 = 0.9;
nu = beta0;

for t = 1:n_iter
    disp(['iteration #' int2str(t)]);
    % update x
    x = support_constraint( y - lambda, mask );
    % update pi0

    % update y
    y = yconstraint( x + pi0, fimg );
    % update lambda
    lambda = lambda + beta0 * (x-y);
end


figure(2001);
imagesc( x );
axis image;
colormap(gray);
colorbar;

figure(2002);
imagesc( y );
axis image;
colormap(gray);
colorbar;

figure(2003);
imagesc( lambda );
colormap(jet);
axis image;
colorbar;
