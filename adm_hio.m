function [x,y, lambda, ers] = adm_hio(fimg, mask, n_iter, img)
rng(911);

Px = @(xin, mask) xin .* (mask & xin > 0);
Py = @(yin, fimg) ifft2( abs(fimg) .* exp( 1j * angle(fft2(yin)) ), 'symmetric' );

z_init = rand(size(fimg));
x_init = Px( z_init, mask );
y_init = Py( z_init, fimg );



%% ADM

x = x_init;
y = y_init;
lambda = zeros(size(x));

beta0 = 0.9;

ers = zeros(2, n_iter);

for t = 1:n_iter
%     disp(['iteration #' int2str(t)]);
    % update x
    x = Px( y - lambda, mask );
    % update y
    y = Py( x + lambda, fimg );
    % update lambda
    lambda = lambda + beta0 * (x-y);
    
    ers(1,t) = er(img, myalign(img,x), mask);
    ers(2,t) = er(img, myalign(img,y), mask);
end

