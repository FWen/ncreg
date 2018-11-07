clear all; clc;close all;

WIDTH = 256;
n = WIDTH^2;   	% signal dimension
m = round(0.4*n);    	% number of measurements

XX(:,:,1) = phantom('Modified Shepp-Logan',WIDTH);
XX(:,:,2) = imresize(double(imread('mri1.jpg')),[WIDTH WIDTH]);

J = randperm(n); J = J(1:m);    % m randomly chosen indices
A = partialDCT(n,m,J); 

figure(1);subplot(121);imagesc(XX(:,:,1));
set(gca,'ytick',[]);set(gca,'xtick',[]);box off;title('Shepp-Logan');
subplot(122);imagesc(XX(:,:,2));
set(gca,'ytick',[]);set(gca,'xtick',[]);box off;title('MRI');

xx =[];
for iX=1:size(XX,3)
    X = XX(:,:,iX);
    xxm1 = xx;

    %Obtain wavelet coeffs
    [C,S] = wavedec2(X, 3, 'haar');
    Norm_C = norm(C);
    x = C'/Norm_C;

    y0 = A*x;
       
    % --noise---------------------------------
    noise = stblrnd(2, 0, 1e-5, 0, m, 1);%Gaussian noise,
        
    y = y0 + noise;
    
    t0 = tic;

    lamdas = logspace(-8, -2, 40);
    % ----L1--------------
    parfor k = 1:length(lamdas)
        [xr, out] = l1_l2_fista(A, y, lamdas(k), x);
        relerr(k) = norm(xr - x)/norm(x);
        xx(:,k)   = xr;
    end
    [~, mi] = min(relerr);
    
    figure(12);semilogy(lamdas,relerr,'r-*');set(gca,'xscale','log');grid;
    
    imgr = waverec2(xx(:,mi)*Norm_C, S, 'haar');
    PSNR(1) = psnr(X, imgr);
    
    figure(2);subplot(2,2,1+(iX-1)*2);
    imagesc(imgr);
    set(gca,'ytick',[]);set(gca,'xtick',[]);box off;
    title(['L1 (', 'PSNR=', num2str(PSNR(1), '%10.2f'), ' dB)']);

    x0 = xx(:,mi);


    % ----Lq--------------
    parfor k = 1:length(lamdas)
        [xr, out] = lq_l2_ista(A, y, lamdas(k), 0.5, x, x0);
        relerr(k) = norm(xr - x)/norm(x);
        xx(:,k)   = xr;
    end
    [~, mi] = min(relerr);
        
    figure(13);semilogy(lamdas,relerr,'r-*');set(gca,'xscale','log');grid;

    imgr = waverec2(xx(:,mi)*Norm_C, S, 'haar');
    PSNR(2) = psnr(X, imgr);
    figure(2);subplot(2,2,2+(iX-1)*2);
    imagesc(imgr);
    set(gca,'ytick',[]);set(gca,'xtick',[]);box off;
    title(['Lq, q=0.5 (', 'PSNR=', num2str(PSNR(2), '%10.2f'), ' dB)']);
    
end

