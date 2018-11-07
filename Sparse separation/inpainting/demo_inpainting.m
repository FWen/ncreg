clear all; clc; close all;

image_files = {'Lenna.tiff'};

HIGHT = 64;
WIDTH = HIGHT;
n = WIDTH*HIGHT;

for img=1:size(image_files,2)

    XX = imresize(imread(cell2mat(image_files(img))),[HIGHT WIDTH]);
    
    A1 = mtxIDCT(n); 
    A2 = mtxEYE(n);

    X = double([reshape(XX(:,:,1),[n,1]), reshape(XX(:,:,2),[n,1]), reshape(XX(:,:,3),[n,1])]);
    x = A1'*X;%true dct coefficients

    y = X;   
    J = randperm(n); 
    J = J(1:round(0.3*n));   %corruption ratio
    y(J,:) = 255;
    y(J(1:round(length(J)/2)),:) = 255;
    y(J(round(length(J)/2)+1:end),:) = 0;

    
    figure(1);subplot(3,3,1+(img-1)*4);imshow(v_2_color_imag(y, HIGHT, WIDTH));
    title(sprintf('Corrupted\n RelErr=%.3f, PSNR=%.2f dB',norm(y-X,'fro')/norm(X,'fro'),psnr(y, X)));
    t0 = tic;

   
    % YALL1--------------------
    n_mu = 15;
    mu     = logspace(-1.5, 1,n_mu);
    relerr_yal = zeros(1,n_mu);
    xx_yal = zeros(n,3,n_mu);   
    parfor k = 1:length(mu); 
        xr = YALL1_admm(A1, y, mu(k), zeros(size(x)), x);
        relerr_yal(k)  = norm(xr-x,'fro')/norm(x,'fro');
        xx_yal(:,:,k) = xr;
    end
    [RelErr(2) mi] = min(relerr_yal);
    x_YALL1(:,:,img) = xx_yal(:,:,mi); 

    figure(1);subplot(3,3,2+(img-1)*4);imshow(v_2_color_imag(idct(x_YALL1(:,:,img)), HIGHT, WIDTH));
    title(sprintf('L1\n RelErr=%.3f, PSNR=%.2f dB',RelErr(2),psnr(idct(x_YALL1(:,:,img)), X)));


    % Lq-Lq-BCD--------------------       
    
    % Initialization. S-ADMM (with q1=q2=1 and \mu=1) and JP solves the same convex formulation, 
    % and have the same accuracy
    x01 = x_YALL1(:,:,img);
    x02 = y - A1*x01;
    
    q1 = 0.7;
    q2 = 0.4;

    t0=tic;
    parfor k = 1:length(mu)
        [xr1] = lq_lq_l2_bcd(A1, A2, y, mu(k), q1, q2, x01, x02, x);
        relerr_bcd(k) = norm(xr1-x,'fro')/norm(x,'fro');
        xx_bcd(:,:,k) = xr1;
    end
    [RelErrs_bcd, mi] = min(relerr_bcd); 
    x_bcd(:,:,img)  = reshape(xx_bcd(:,:,mi),[n*3,1]); 

    figure(1);subplot(1,3,3+(img-1)*4);
    imshow(v_2_color_imag(idct(reshape(x_bcd(:,:,img),[n,3])), HIGHT, WIDTH));
    title(sprintf('Lq (q1=%.1f, q2=%.1f)\n RelErr=%.3f, PSNR=%.2f dB',q1,q2,RelErrs_bcd,psnr(idct(reshape(x_bcd(:,:,img),[n,3])), X)));

end
