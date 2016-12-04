load('lena_101.mat');

img = im2double( img_padded );
fimg = abs( fft2( img ) );

% figure(1001);
% imshow(img);

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

sig = 1;

G = exp( -(RR).^2/(2*sig) );

% figure(3001);
% imagesc(G);
% axis image;
% colormap(jet);
% colorbar;


img_lrs = ifft2( fft2(img) .* ifftshift(G), 'symmetric');

figure(3002);
imshow(img_lrs);




%% Initialization

rng(5566);

z_init = rand(size(img));

x_init = Px( z_init, mask );
y_init = Py( z_init, fimg );

% figure(1011);
% imshow(x_init);
% 
% figure(1012);
% imshow(y_init);

%% ADMM pre-conputed functions

max_pcg_iter = 100;
lambda = 1;
rho = 25;

dx = [0 0 0; 0 -1 1; 0 0 0];
dy = dx.';


kappa = lambda / rho;
p2o = @(h) psf2otf(h, size(img));

% D function family
dxFT = p2o(dx);
dxTFT = conj(p2o(dx)); 
dyFT = p2o(dy);
dyTFT = conj(p2o(dy));

Dfun    = @(x) cat(3, ifft2( fft2(x) .* dxFT), ifft2( fft2(x) .* dyFT ) );
Dtfun   = @(x) ifft2( fft2(x(:,:,1)) .* dxTFT + fft2(x(:,:,2)) .* dyTFT ); 
DtDfun  = @(x) ifft2 ( fft2 (x) .* dxFT.* dxTFT + fft2 (x) .* dyFT.* dyTFT ); 


% A function family
Afun  = @(x) Py(x, fimg);
Atfun = @(x) conj(Py(x, (fimg)));
% Atfun = @(x) x;

A_tilde = @(x) reshape(...
    Atfun(Afun(reshape(x, size(img)))) + rho * DtDfun(reshape(x, size(img))),...
    [numel(img) 1]);



%% ADM

n_iter = 5;

x = x_init;
y = y_init;

z = rand(size(x,1), size(x,2), 2);
u = rand(size(x,1), size(x,2), 2);

lambda = zeros(size(x));

beta0 = 0.8;
nu0 = 0.75;
ers = zeros(1, n_iter);

% img_lrs = 1 * mask;

b = img;%_lrs;

for t = 1:n_iter
    disp(['iteration #' int2str(t)]);
    % update x
    v = z - u;
    b_tilde = reshape(...
        Atfun(b) + rho * Dtfun(v),...
        [numel(img) 1]);
%     x = pcg(A_tilde, b_tilde, 1e-12, max_pcg_iter);
%     x = reshape(x, size(img));
%     x = Px(x, mask);
    x1 = fft2( Py(b, conj(fimg))) + rho * (dxTFT.*fft2(v(:,:,1)+dyTFT.*fft2(v(:,:,2)) ) ) ; 
    x2 = 0 + rho * (dxTFT.*dxFT+dyTFT.*dyFT );
    x = ifft2( x1 ./ x2, 'symmetric' );
    % update z
    v = Dfun(x) + u;
    z( v > kappa ) = v( v > kappa ) - kappa;
    z( v < -kappa ) = v( v < -kappa ) + kappa;
    z( abs(v) <= kappa ) = 0;
    z( cat(3, ~mask, ~mask) ) = 0;
    % update u
    u = u + Dfun(x) - z;
    
    u( cat(3, ~mask, ~mask) ) = 0;
    
    
    
    ers(1,t) = er(img, reshape(x, size(img)), mask);
end

%% plot

er_x = er(img, x, mask);

x = reshape(x, size(img));

figure(2001);
imagesc( x );
axis image;
colormap(gray);
colorbar;
title(['E_R = ' num2str(er_x)]);

% figure(2005);
% plot(1:n_iter, ers, 'LineWidth', 2);
% set(gca, 'FontSize', 18);
% legend('x');
figure(2006);
semilogy(1:n_iter, ers, 'LineWidth', 2);
set(gca, 'FontSize', 18);
legend('x');


