function [x,out] = YALL1_admm(A,y,lamda,x0,xtrue);
% YALL1_admm solves
%
%   minimize || Ax - y ||_1 + \lamda || x ||_1
%
% Inputs:
%	A: sensing matrix
%	y: CS data
%	lamda: regularization parameter 
%	x0: initialization 
%	xtrue: for debug, for calculation of errors
% Outputs
%	x: the recovery
%	out.e: the error with respect to the true
%	out.et: time index


%Convergence setup
max_iter = 300;
ABSTOL   = 1e-6;

[m, L] = size(y);

if(isobject(A))
    n = A.n;
else
    n=size(A,2);
end

rho_cov = 1e2;
rho  = 1e-2;
    
%Initialize
if nargin<4
	x = zeros(n,L);
else
	x = x0;
end;

isquiet = 0;
if nargin<5
    isquiet = 1;
end

w = zeros(m,L); 

out.e  = [];
out.et = [];
tic;

for i = 1 : max_iter
	
    if (rho<rho_cov)
        rho = rho*1.02;
    end
    
    xm1 = x;
    
	%v-step
	tv = A*x-y-w/rho;
    v = repmat(prox_Lq_m(tv,1,lamda*rho),[1,L]) .* tv;
    
    %x-step
    tao = 0.99; %for orthonormal A
    z = x - tao*(A'*(A*x - y - v - w/rho)); 
    x = repmat(prox_Lq_m(z,1,rho/tao),[1,L]) .* z;
    
	%u-step
    Ax = A*x;
	w = w - rho*(Ax - y - v);

    if ~isquiet
        out.e  = [out.e norm(x-xtrue)/norm(xtrue)];
        out.et = [out.et toc];
    end
    
    %terminate when residuals are small
    if (norm(x-xm1)< sqrt(n)*ABSTOL && norm(Ax-y-v)< sqrt(n)*ABSTOL) 
        break;
    end
end

end
