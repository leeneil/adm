% 2-D HIO written by Po-Nan Li @ Academia Sinica 2012
function [R, efs, ers] = hio2d_lr_rs(Fabs, S, n, img_lr, img_hr, varargin) % Fabs, S, n, unknown, alpha
    
    rng(911);

    if isempty(varargin)
        unknown = false(size(Fabs));
    else
        unknown = varargin{1};
    end
    % OSS module
    if length(varargin) > 1 && ~isempty( varargin{2} )
        alpha = varargin{2};
        disp(['OSS is on. alpha = ' num2str(alpha)]);
        oss = true;
        x = -round((length(Fabs)-1)/2):round((length(Fabs)-1)/2);
        [X, Y] = meshgrid(x, x);
        W = exp(-0.5 .* (X./alpha).^2) .* exp(-0.5 .* (Y./alpha).^2);
        W = ifftshift(W);
    else
        oss = false;
    end
    % solve unknown pixels in data
    
    
    beta1 = 0.9;
    gamma = 1;
    
    % convergence parameters
    cvg_cutoff1 = 0.0001;
    cvg_cutoff2 = 0.0005;
    cvg_stop = 10;
    cvg_count = 0;
    cvg_s1_img = zeros( size( Fabs ) );
    
    
    % generate random initial phases
    if sum(imag(Fabs(:))) == 0
        ph_init = rand(size(Fabs));
        ph_init = angle(fft2(ph_init));
        F = Fabs .* exp(1j.*ph_init);
    else
        F = Fabs;
    end
    
    F0 = abs(F); 
    previous = ifft2(F, 'symmetric');
    
    
    efs = ones(2, n);
    ers = ones(2, n);
    
    % ================ iterations ==================================
    for t = 1:n
        if mod(t-1, 100) == 0 && n >= 500
            disp(['step ' int2str(t)]);
        end
        rs = ifft2(F, 'symmetric'); % real space version
        cond1 = ~S | (rs<0);
        
        rs(cond1) = previous(cond1) - beta1 .* rs(cond1);
        
        ers(1,t) = er(img_hr, myalign(img_hr, rs), S);
        
%         hr_er = er(img_hr, rs, S);
%         ers(t) = hr_er;

        rs_ef = rs;
        rs_ef(cond1) = 0;
        efs(1,t) = ef(F0, fft2(rs_ef), unknown);
        
        
        cm = img_lr ./ rs;
        cm(isnan(cm)) = 0;
        if t > 1
%             rs = rs + hr_er *  (rs .* cm - rs);
            rs = rs + gamma * efs(t) *  (rs .* cm - rs);
        end
        
        ers(2,t) = er(img_hr, myalign(img_hr, rs), S);
            
        cvg_er = er(previous, rs, S);
        if cvg_er < cvg_cutoff1
            if cvg_count == 0
                cvg_s1_img = rs;
            end
            cvg_count = cvg_count + 1;
        else
            cvg_count = 0;
        end
        
        previous = rs;
        if oss
            rs_oss = ifft2(fft2(rs) .* W, 'symmetric');
            rs(~S) = rs_oss(~S);
        end
        F2 = fft2(rs); % .* exp(-1j.*(U+V));
        
        rs_ef = rs;
        rs_ef(cond1) = 0;
        efs(2,t) = ef(F0, fft2(rs_ef), unknown);
        
        F = F0 .* exp(1j.*angle(F2));
        F(unknown) = F2(unknown);
        if cvg_count >= cvg_stop
            if er(cvg_s1_img, rs, S) < cvg_cutoff2
                disp(['Early termination at iteration #' int2str(t) '.']);
                disp(['Change = ' num2str(cvg_er)]);
                break;
            else
                cvg_count = 0;
            end
        end
    end
        % ================ iterations ends here  ==================================
    R = ifft2(F, 'symmetric');
end