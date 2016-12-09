n_iter = 100;

wn = 0.0001;

% load('lena_101.mat');
load('mimivirus_img_105.mat'); % mimivirus

img = im2double(img_padded);
fimg = abs(fft2(img));

fimg = sqrt ( wn^-1 * poissrnd(wn*fimg.^2 ) );


timers = zeros(5,1);

figure(3001);
imagesc( fftshift( log10(abs(fimg).^2) ) );
axis image;
axis off;
colormap(jet);
colorbar;


%% generate LR image

xx = linspace(-0.5, 0.5, size(img, 2));
yy = linspace(-0.5, 0.5, size(img, 1));

[XX, YY] = meshgrid(xx, yy);
RR = sqrt( XX.^2 + YY.^2 );

sig = 0.01;
G = exp( -(RR).^2/(2*sig) );

% figure(3001);
% imagesc(G);
% axis image;
% colormap(jet);
% colorbar;


img_lrs = ifft2( fft2(img) .* ifftshift(sqrt(G)), 'symmetric');

figure(3002);
imshow(img_lrs);
imwrite(img_lrs, ['minivirus_lrs_sig' int0str(10000*sig, 5) '.png']);

%% test 1: ADM-HIO

tic;
[x1,y1,lambda1,ers1] = adm_hio(fimg, mask, n_iter, img, fimg==0);
toc
timesr(1) = toc;

figure(3011);
plot(1:n_iter, ers1, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x', 'y');

%% test 2: ADM-HIO-LR

tic;
[x2,y2,lambda2,ers2] = adm_hio_lr(fimg, mask, img_lrs, n_iter, img, fimg==0);
toc
timers(2) = toc;

figure(3012);
plot(1:n_iter, ers2, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x', 'y');

%% test 3: ADM-HIO-LR-TV

%% test 4: HIO-LR-RS

tic;
[x4, efs4, ers4] = hio2d_lr_rs(fimg, mask, n_iter, img_lrs, img, fimg==0);
toc
times(4) = toc;

figure(3014);
plot(1:n_iter, ers4, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x1', 'x2');

%% test 5: HIO

tic;
[x5, efs5, ers5] = hio2d(fimg, mask, n_iter, fimg==0, [], img);
toc
timers(5) = toc;

figure(3015);
plot(1:n_iter, ers5, 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('x');
%% comparison

t = 1:n_iter;

figure(5001);
plot(t, ers5, t, ers1(1,:), t, ers4(1,:), t, ers2(1,:), 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('HIO','ADM-HIO', 'HIO-LR', 'ADM-HIO-LR');

figure(5002);
semilogy(t, ers5, t, ers1(1,:), t, ers4(1,:), t, ers2(1,:), 'LineWidth', 2);
set(gca, 'FontSize', 16);
xlabel('iteration number');
ylabel('E_R');
legend('HIO (Fienup)','ADM-HIO (Marchesini)', 'HIO-LR (Li)', 'ADM-HIO-LR (Li & Huang)');

save(['mimivirus_lrs_sig' int0str(10000*sig, 5) '_nos' int0str(10000*wn, 7) '.mat'], 'x*', 'y*', 'ers*', 'efs*','timers', 't');
