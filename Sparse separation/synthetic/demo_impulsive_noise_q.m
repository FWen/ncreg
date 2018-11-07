clear all; clc;

N = 256;
M = 100;

A1 = dctmtx(N);
J = randperm(N); 
A1 = A1(J(1:M),:);

A2 = eye(M);

K = 20;

x1 = 5*SparseVector(N,K);
x2 = stblrnd(1, 0, 1e-3, 0, M, 1);

y  = A1*x1 + x2;
figure(1);plot(x2);ylabel('Amplitude');xlabel('Sample index');xlim([1,length(x2)])
        
mu = logspace(-2, 1,15);


% S-ADMM with q1=q2=1 and \mu=1, for initialization
[x01,x02,~] = lq_lq_admm(A1, A2, y, 1, 1, 1, zeros(N,1), zeros(M,1), x1);


% Lq-Lq-BCD and Lq-Lq-ADMM--------------------
qs = 0:0.1:1;

for l1=1:length(qs)
    for l2=1:length(qs)
        parfor k = 1:length(mu)
           
            [xr1,xr2,out1] = lq_lq_l2_bcd(A1, A2, y, mu(k), qs(l1), qs(l2), x01, x02, x1);
            relerr_bcd(k) = norm(xr1-x1)/norm(x1);
            
            [xr1,xr2,out2] = lq_lq_l2_admm(A1, A2, y, mu(k), qs(l1), qs(l2), x01, x02, x1);
            relerr_admm(k) = norm(xr1-x1)/norm(x1);

        end
        
        [RelErrs_bcd(l1,l2), mi] = min(relerr_bcd);       
        [RelErrs_admm(l1,l2), mi] = min(relerr_admm);        

        sprintf('q1=%.1f and q2=%.1f completed',qs(l1),qs(l2))
    end
end

figure(2);
subplot(121);v0=-55;v1=-20;
[w1, e1] = min(RelErrs_bcd);[~, lo] = min(w1); ko = e1(lo);
contourf(qs,qs,20*log10(RelErrs_bcd));  title('BCD');
colorbar; xlabel('q_2');ylabel('q_1');hold on;%set(gca, 'CLim', [v0, v1]);

subplot(122);
[w1, e1] = min(RelErrs_admm);[~, lo] = min(w1); ko = e1(lo);
contourf(qs,qs,20*log10(RelErrs_admm));  title('ADMM');
colorbar; xlabel('q_2');ylabel('q_1');hold on;%set(gca, 'CLim', [v0, v1]);
