clear all;clc;close all;

N = 100;
M = 40;

A = randn(M,N);
A = orth(A')'; 

MC = 300;
Ks = [2, 5, 10];
SNRs = 0:5:40;

for n=1:length(Ks)
    for snr=1:length(SNRs)

        for l = 1:MC
            disp(['Sparsity: ', num2str(Ks(n)), ',   SNR: ', num2str(SNRs(snr)), ',   Time: ', num2str(l)]);

            t0 = tic;
            x  = SparseVector(N, Ks(n));
            y0 = A*x;

            % --Gaussian noise-----------------------------
            SNR   = SNRs(snr);  % dB
            noise = randn(M,1);
            noise = noise/std(noise) *10^(-SNR/20)*std(y0);  
            y = y0 + noise;

            lamdas = logspace(-5.5,1,20);

            %---- L1 ------------------------------------
            for k = 1:length(lamdas)
                [x_rec1, out1] = lq_l2_fista(A, y, lamdas(k), 1, x);
                relerr(1,k)   = norm(x_rec1 - x)/norm(x);
                xx(:,k)       = x_rec1;
            end
            [mv mi] = min(relerr(1,:));
            x_0 = xx(:,mi); 
           

            %---- L0, Lq, SCAD -------------------------
            for k = 1:length(lamdas)
                % LqLS-FISTA (q=0.7)
                [x_rec2, out2] = lq_l2_fista(A, y, lamdas(k), 0.7, x, x_0);
                relerr(2,k)   = norm(x_rec2 - x)/norm(x);

                % LqLS-FISTA (q=0.5)
                [x_rec3, out3] = lq_l2_fista(A, y, lamdas(k), 0.5, x, x_0);
                relerr(3,k)    = norm(x_rec3 - x)/norm(x);
                
                % LqLS-FISTA (q=0.2)
                [x_rec4, out4] = lq_l2_fista(A, y, lamdas(k), 0.2, x, x_0);
                relerr(4,k)   = norm(x_rec4 - x)/norm(x);               

                % LqLS-FISTA (q=0.0)
                [x_rec5, out5] = lq_l2_fista(A, y, lamdas(k), 0.0, x, x_0);
                relerr(5,k)    = norm(x_rec5 - x)/norm(x);

                % SCAD-ISTA
                [x_rec6, out6] = scad_ista(A, y, lamdas(k), x, x_0);
                relerr(6,k)    = norm(x_rec6 - x)/norm(x);
            end
%           figure(2);
%           semilogy(lamdas,relerr(1,:),'r--*',lamdas,relerr(2,:),'b--*',lamdas,relerr(3,:),'g--*',...
%                      lamdas,relerr(4,:),'k--*',lamdas,relerr(5,:),'c--*',lamdas,relerr(6,:),'c--+','linewidth',1); set(gca,'xscale','log');
            
            RelErr(l,:) = min(relerr')';

            toc(t0)
        end

        aver_Err(snr,(1:6)+(n-1)*6) = mean(RelErr);
    end
   
    figure(1); subplot(1,3,n);
    semilogy(SNRs,aver_Err(:,1+(n-1)*6),'r-',SNRs,aver_Err(:,2+(n-1)*6),'b--',SNRs,aver_Err(:,3+(n-1)*6),'b-.',...
             SNRs,aver_Err(:,4+(n-1)*6),'b:+',SNRs,aver_Err(:,5+(n-1)*6),'g--',SNRs,aver_Err(:,6+(n-1)*6),'k:*');
    legend('L1','Lq (q=0.7)','Lq (q=0.5)','Lq (q=0.2)','Hard','SCAD','Location','Best');grid off;
    ylabel('Averaged relative error of recovery'); xlabel('SNR (dB)');
        
end

