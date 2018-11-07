clear all; clc; close all;

WIDTH = 500;
HIGHT = 318;
n = WIDTH*HIGHT;   	% signal dimension

X0 = imresize(imread('turtle.png'),[HIGHT WIDTH]);
X1 = imresize(imread('masked_transport.png'),[HIGHT WIDTH]);%corrupted by text

figure(1);subplot(2,2,1);imshow(X0);set(gcf,'outerposition',get(0,'screensize'));title('Original');

A1 = mtxIDCT(n); 
A2 = mtxEYE(n);

X = double([reshape(X0(:,:,1),[n,1]), reshape(X0(:,:,2),[n,1]), reshape(X0(:,:,3),[n,1])]);
x = A1'*X; %true dct coefficients
y = double([reshape(X1(:,:,1),[n,1]), reshape(X1(:,:,2),[n,1]), reshape(X1(:,:,3),[n,1])]); 

figure(1);subplot(2,2,2);imshow(v_2_color_imag(y, HIGHT, WIDTH));
title(sprintf('Corrupted\n RelErr=%.3f, PSNR=%.2f dB',norm(y-X,'fro')/norm(X,'fro'),psnr(y, X)));

t0=tic;

% YALL1--------------------
mu_min = 1e-1; mu_max = 1e1; n_mu = 15;
mu     = logspace(log10(mu_min), log10(mu_max),n_mu);
relerr_yal = zeros(1,n_mu);
xx_yal = zeros(n,3,n_mu);   
parfor k = 1:length(mu); 
    xr = YALL1_admm(A1, y, mu(k), zeros(size(x)), x);
    relerr_yal(k)  = norm(xr-x,'fro')/norm(x,'fro');
    xx_yal(:,:,k) = xr;
end
[RelErr(2) mi] = min(relerr_yal);
x_YALL1 = xx_yal(:,:,mi); 

figure(1);subplot(2,2,3);imshow(v_2_color_imag(idct(x_YALL1), HIGHT, WIDTH));
title(sprintf('YALL1\n RelErr=%.3f, PSNR=%.2f dB',RelErr(2),psnr(idct(x_YALL1), X)));

sprintf('YALL1 completed, elapsed time: %.1f seconds',toc(t0))
t0=tic;


% Lq-Lq-BCD and Lq-Lq-ADMM--------------------
x01 = x_YALL1;
x02 = y - A1*x_YALL1;

parfor k = 1:length(mu)
    [xr1] = lq_lq_l2_bcd(A1, A2, y, mu(k), 0.7, 0.4, x01, x02, x);
    relerr_bcd(k) = norm(xr1-x,'fro')/norm(x,'fro');
    xx_bcd(:,:,k) = xr1;
end
[RelErrs_bcd, mi] = min(relerr_bcd); 
x_bcd  = reshape(xx_bcd(:,:,mi),[n*3,1]); 
PSNR_bcd = psnr(idct(reshape(x_bcd,[n,3])), X);

figure(1);subplot(2,2,4);
imshow(v_2_color_imag(idct(reshape(x_bcd,[n,3])), HIGHT, WIDTH));
title(sprintf('BCD (q1=0.7, q2=0.4)\n RelErr=%.3f, PSNR=%.2f dB',RelErrs_bcd,PSNR_bcd));
