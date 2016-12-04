% 2-D HIO written by Po-Nan Li @ Academia Sinica 2012
function [R, efs, ers] = hio2d_lr_fs(Fabs, S, n, unknown, img_hs) % Fabs, S, n, unknown, alpha
    
    rng(911);
    
    if isempty(unknown)
        unknown = false(size(Fabs));
    end
    
    
    
    beta1 = 0.9;
    
    % convergence parameters
    cvg_cutoff1 = 0.0001;
    cvg_cutoff2 = 0.00050;
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
    ers = zeros(n, 1);
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
        
        F2 = fft2(rs); % .* exp(-1j.*(U+V));
        rs_er = rs;
        rs_er(cond1) = 0;
        
        efs(t) = ef(F0, fft2(rs_er), unknown);
        ers(t) = er(img_hs, myalign(img_hs, rs_er), S);
        F = F0 .* exp(1j.*angle(F2));
        F(unknown) = efs(t) * F(unknown) + (1-efs(t)) * F2(unknown);
%         F(unknown) = ( efs(t) * abs(F(unknown)) + (1-efs(t)) * abs(F2(unknown)) ) .* exp(1j.*angle(F2(unknown)));
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