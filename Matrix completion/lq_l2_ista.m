function [X, out] = lq_l2_ista(q, M, P, lamda, Xtrue, X0);

[m,n] = size(M);

% Lipschitz constant
L = 1.1;

%Convergence setup
MAX_ITER = 2000;
ABSTOL = 1e-5;

%Initialize
if nargin<6
	X = zeros(m,n);
else
    X = X0;
end

out.et = [];
tic;

for iter = 1 : MAX_ITER

    Xm1 = X;	
    V = X - (1/L)*((X - M).*P);

    [S,V,D] = svd(V);
    
    v = shrinkage_Lq(diag(V), q, lamda, L);
    
    X = S*diag(v)*D';

    
%     out.e  = [out.e, norm(X-Xtrue,'fro')/norm(Xtrue,'fro')];
        
    %Check for convergence
    if norm(X-Xm1,'fro')<ABSTOL*sqrt(m*n)
        break;
    end

end


