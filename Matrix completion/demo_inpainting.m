clear all; clc;  close all;

WIDTH = 512;
X = imresize(double(imread('boat.png')),[WIDTH, WIDTH]);
[m,n] = size(X);

[S,V,D] = svd(X);
v = (diag(V));
% figure(2);
% plot(1:n,v/max(v),'-'); xlim([1,n]);

figure(3);subplot(121);
semilogy(1:m,v,'r-','linewidth',1); xlim([1,512]);%ylim([0.1,1e5]);
ylabel('Singular value');title('Not strictly low-rank');
legend('NSLR');


v(round(WIDTH*0.15):end) = 0;
X = S*diag(v)*D';
% figure(1);subplot(2,2,2);
% imshow(uint8(X2));

subplot(122);
semilogy(1:m,v,'b-','linewidth',2); xlim([1,512]);%ylim([0.1,1e5]);
% ylabel('RelErr');xlabel('q');
legend('SLR');title('Strictly low-rank')


J = randperm(m*n); 
J = J(1:round(0.5*m*n));    % sampling ratio
P = ones(m*n,1);
P(J) = 0;
P = reshape(P,[m,n]);   % Projection Matrix

SNR = 40;
M = X(:);
noise = randn(m*n,1);
noise = noise/std(noise) *10^(-SNR/20)*std(M);  
M = M + noise;
M = reshape(M,[m,n]).*P; % partial observation

figure(1); subplot(1,4,1);
imshow(uint8(M));
title('Partial observation (50%)');


% L1  --------------------
lamdas = logspace(-1, 3,20);
parfor k = 1:length(lamdas); 
    [Xr, out]   = lq_l2_ista(1, M, P, lamdas(k), X, zeros(size(M)));
    relerr1(k)  = norm(Xr-X,'fro')/norm(X,'fro');
    X_l1(:,:,k) = Xr;
end
[RelErr1, mi] = min(relerr1);
X_L1 = X_l1(:,:,mi); 

subplot(1,4,2);
imshow(uint8(X_L1));
title(sprintf('Soft (RelErr=%.5f)', RelErr1));

% L0  --------------------
parfor k = 1:length(lamdas); 
    [Xr, out2]  = lq_l2_ista(0, M, P, lamdas(k), X, X_L1);
    relerr3(k)  = norm(Xr-X,'fro')/norm(X,'fro');
    X_l0(:,:,k) = Xr;
end

[RelErr3, mi] = min(relerr3);
X_L0 = X_l0(:,:,mi);
subplot(1,4,4);
imshow(uint8(X_L0));
title(sprintf('Hard (RelErr=%.5f)', RelErr3));

% Lq  --------------------
qs = 0.1:0.1:0.9;
for kq=1:length(qs)
    qs(kq)
    parfor k = 1:length(lamdas)       
        [Xr,~]     = lq_l2_ista(qs(kq), M, P, lamdas(k), X, X_L1);
        relerr(k)  = norm(Xr-X,'fro')/norm(X,'fro');
        X_lq(:,:,k) = Xr;
    end
    [RelErr, mi]  = min(relerr);
    imgrs(:,:,kq) = X_lq(:,:,mi);
    RelErrs(kq)   = RelErr;
end

% figure(4);semilogy(lamdas,relerr1,'r-o',lamdas,relerr3,'g-x');
% set(gca,'xscale','log');

[min_RelErr, mi] = min(RelErrs);

subplot(1,4,3);
imshow(uint8(imgrs(:,:,mi)));
title(['Lq (best q=', num2str(qs(mi), '%.1f'), ', ', 'RelErr=', num2str(min_RelErr, '%.5f'), ')']);

figure;
semilogy(0:0.1:1,[RelErr3, RelErrs, RelErr1],'--*');
xlabel('q');ylabel('Relative error of recovery');

