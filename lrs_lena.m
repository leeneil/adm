load('lena_101.mat');

img = im2double( img_padded );
fimg = abs( fft2( img ) );

figure(1001);
imshow(img);

% figure(1002);
% imagesc( log10( fftshift(fimg.^2) ) );
% axis image;
% colormap(jet);
% colorbar;

Px = @(xin, mask) xin .* (mask & xin > 0);
Py = @(yin, fimg) ifft2( abs(fimg) .* exp( 1j * angle(fft2(yin)) ), 'symmetric' );

%% generate low-res

xx = linspace(-0.5, 0.5, size(img, 2));
yy = linspace(-0.5, 0.5, size(img, 1));

[XX, YY] = meshgrid(xx, yy);
RR = sqrt( XX.^2 + YY.^2 );

sig = 0.01;

G = exp( -(RR).^2/(2*sig) );

figure(3001);
imagesc(G);
axis image;
colormap(jet);
colorbar;


img_lrs = ifft2( fft2(img) .* ifftshift(G), 'symmetric');

figure(3002);
imshow(img_lrs);




%% Initialization

rng(5566);

% z_init = rand(size(img));
z_init = zeros(size(img));

x_init = support_constraint( z_init, mask );
y_init = Py( z_init, fimg );

% figure(1011);
% imshow(x_init);
% 
% figure(1012);
% imshow(y_init);

%% ADM

n_iter = 50;

x = x_init;
y = y_init;
lambda = zeros(size(x));

beta0 = 0.8;
nu0 = 0.75;
ers = zeros(2, n_iter);

% img_lrs = 1 * mask;


for t = 1:n_iter
    disp(['iteration #' int2str(t)]);
    % update x
%     x = support_constraint( img_lrs - lambda, mask );
%     x = Px( img_lrs - lambda, mask );
    x = img_lrs - lambda;
    % update y
    y = Py( x + lambda, fimg );
    % update lambda
    lambda = lambda + beta0 * (x-y);
    
    ers(1,t) = er(img, x, mask);
    ers(2,t) = er(img, y, mask);
end

%% plot

er_x = er(img, x, mask);
er_y = er(img, y, mask);

figure(2001);
imagesc( x );
axis image;
colormap(gray);
colorbar;
title(['E_R = ' num2str(er_x)]);

figure(2005);
plot(1:n_iter, ers, 'LineWidth', 2);
set(gca, 'FontSize', 18);
legend('x', 'y');
figure(2006);
semilogy(1:n_iter, ers, 'LineWidth', 2);
set(gca, 'FontSize', 18);
legend('x', 'y');


figure(2002);
imagesc( y );
axis image;
colormap(gray);
colorbar;
title(['E_R = ' num2str(er_y)]);

figure(2003);
imagesc( lambda );
colormap(jet);
axis image;
colorbar;