% 2-D HIO written by Po-Nan Li @ Academia Sinica 2012
function [R, efs, ers, cvg_ers, F0, F2] = hio2d(Fabs, S, n, varargin) % Fabs, S, n, unknown, alpha
    
    rng(911);

    if isempty(varargin)
        unknown = false(size(Fabs));
    else
        unknown = varargin{1};
    end
    
    otf_er = false;
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
    
    if length(varargin) > 2 && ~isempty( varargin{3} )
        otf_er = true;
        ers = zeros(1, n);
        img_tmp = varargin{3};
    end
    
    % solve unknown pixels in data
    
    
    beta1 = 0.9;
    
    % convergence parameters
    cvg_cutoff1 = 0.000001;
    cvg_cutoff2 = 0.0000050;
    cvg_stop = 20;
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
    
    efs = zeros(n, 1);
    cvg_ers = zeros(n, 1);
    
    F0 = abs(F); 
    previous = ifft2(F, 'symmetric');
    
    % ================ iterations ==================================
    for t = 1:n
        if mod(t-1, 100) == 0 && n >= 500
            disp(['step ' int2str(t)]);
        end
        rs = ifft2(F, 'symmetric'); % real space version
        cond1 = ~S | (rs<0);
        rs(cond1) = previous(cond1) - beta1 .* rs(cond1);
        cvg_er = er(previous, rs, S);
        cvg_ers(t) = cvg_er;
        if cvg_er < cvg_cutoff1
            if cvg_count == 0
                cvg_s1_img = rs;
            end
            cvg_count = cvg_count + 1;
        else
            cvg_count = 0;
        end
        
        previous = rs;
        
        if otf_er
            ers(t) = er(img_tmp, myalign(img_tmp, rs), S);
        end
        
        if oss
            rs_oss = ifft2(fft2(rs) .* W, 'symmetric');
            rs(~S) = rs_oss(~S);
        end
        F2 = fft2(rs); % .* exp(-1j.*(U+V));
        rs_er = rs;
        rs_er(cond1) = 0;
        efs(t) = ef(F0, fft2(rs_er), unknown);
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
%     
    cond1 = ~S | (R<0);
    R(cond1) = 0;
    ef00 = ef(F0, fft2(R), unknown)
end