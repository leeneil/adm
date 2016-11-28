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

x_init = rand(size(img));
y_init = rand(size(img));


x_init( ~mask ) = 0;

figure(1011);
imshow(x_init);

y_init = ifft2( abs(fimg) .* exp( 1j * angle(fft2(y_init)) ), 'symmetric' );

figure(1012);
imshow(y_init);

n_iter = 500;

x = x_init;
y = y_init;
lambda = zeros(size(x));
pi0 = zeros(size(x));

beta0 = 1;
nu = beta0;

for t = 1:n_iter
    disp(['iteration #' int2str(t)]);
    % update x
    x = y - lambda;
    x( ~mask ) = 0;
    x( x < 0 ) = 0;
    % update pi0
    
    % update y
    y = x + pi0;
    y = ifft2( abs(fimg) .* exp( 1j * angle(fft2(y)) ), 'symmetric' );
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

