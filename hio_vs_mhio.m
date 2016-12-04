n_iter = 5000;

load('lena_101.mat');

img = im2double(img_padded);
fimg0 = abs(fft2(img));


figure(3001);
imagesc( fftshift( log10(abs(fimg0).^2) ) );
axis image;
axis off;
colormap(jet);
colorbar;

%% simulate BS

bs = 25;

k = (-(size(img,1)-1)/2):((size(img,1)-1)/2);
pc = (size(img,1)+1)/2;

fimg_bs = fftshift(fimg0);
fimg_bs( (pc-bs):(pc+bs), (pc-bs):(pc+bs) ) = 0;
fimg_bs = ifftshift( fimg_bs );


%% generate LR image

xx = linspace(-0.5, 0.5, size(img, 2));
yy = linspace(-0.5, 0.5, size(img, 1));

[XX, YY] = meshgrid(xx, yy);
RR = sqrt( XX.^2 + YY.^2 );

sig = 0.001;
G = exp( -(RR).^2/(2*sig) );

% figure(3001);
% imagesc(G);
% axis image;
% colormap(jet);
% colorbar;

fimg_lr = fftshift( fft2(img) .* sqrt(ifftshift(G)) );

img_lrs = ifft2( ifftshift( fimg_lr ), 'symmetric');

figure(3002);
imshow(img_lrs);
imwrite(img_lrs, ['lena_lrs_sig' int0str(10000*sig, 5) '.png']);

fimg_bs_lr = fftshift(fimg0);
fimg_bs_lr( (pc-bs):(pc+bs), (pc-bs):(pc+bs) ) = fimg_lr( (pc-bs):(pc+bs), (pc-bs):(pc+bs) );
fimg_bs_lr = ifftshift( fimg_bs_lr );
fimg_lr = ifftshift( fimg_lr );

%% test 1: ADM-HIO

% tic;
% [x1,y1,lambda1,ers1] = adm_hio(fimg, mask, n_iter, img);
% toc
% times(1) = toc;
% 
% figure(3011);
% plot(1:n_iter, ers1, 'LineWidth', 2);
% set(gca, 'FontSize', 16);
% xlabel('iteration number');
% ylabel('E_R');
% legend('x', 'y');

%% test 2: ADM-HIO-LR

% tic;
% [x2,y2,lambda2,ers2] = adm_hio_lr(fimg, mask, img_lrs, n_iter, img);
% toc
% times(2) = toc;
% 
% figure(3012);
% plot(1:n_iter, ers2, 'LineWidth', 2);
% set(gca, 'FontSize', 16);
% xlabel('iteration number');
% ylabel('E_R');
% legend('x', 'y');

%% test 3: ADM-HIO-LR-TV

%% test 4: HIO-LR-RS

tic;
[x4, efs4, ers4] = hio2d_lr_rs(fimg_bs, mask, n_iter, img_lrs, img, fimg_bs==0);
toc


figure(3014);
plot(1:n_iter, ers4, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x1', 'x2');

%% test 5: HIO-LR-FS

tic;
[x5, efs5, ers5] = hio2d_lr_fs(fimg_bs_lr, mask, n_iter, fimg_bs==0, img);
toc


figure(3015);
plot(1:n_iter, ers5, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x');

%% test 6: HIO

tic;
[x6, efs6, ers6] = hio2d(fimg_bs, mask, n_iter, fimg_bs==0, [], img);
toc

figure(3016);
plot(1:n_iter, ers6, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x');

%% test 7: HIO+LR

tic;
[x7, efs7, ers7] = hio2d(fimg_bs_lr, mask, n_iter, fimg_bs==0, [], img);
toc

figure(3017);
plot(1:n_iter, ers7, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x');
%% comparison

t = 1:n_iter;

figure(5001);
plot(t, ers6, 'k', t, ers7, 'k--', t, ers4(1,:), 'r', t, ers5(:), 'b', 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('HIO', 'HIO+', 'HIO-LR-RS', 'HIO-LR-FS');

figure(5002);
semilogy(t, ers6, 'k', t, ers7, 'k--', t, ers4(1,:), 'r', t, ers5(:), 'b', 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('HIO', 'HIO+', 'HIO-LR', 'HIO-LR-FS');

save(['lena_hio_sig' int0str(10000*sig, 5) '_bs' int0str(bs,3) '.mat'], 'x*', 'y*', 'ers*', 'efs*', 't');
