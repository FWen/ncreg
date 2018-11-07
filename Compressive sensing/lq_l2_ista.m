function [x,out] = lq_l2_ista(A,y,lamda,q,xtrue,x0);
% LqLS_ista solves
%
%   minimize || Ax - y ||_2^2 + \lambda || x ||_q^q
%
% Inputs
%	A,y,lambda: CS variables
%	0<=q<=1
%	xtrue: for debug, for calculation of errors
%   x0: intialization
% Outputs
%	x: the CS recovery
%	out.e: the error with respect to the true
%	out.et: time index

Aty = A'*y;
if(isobject(A))
    m=A.m;
    n=A.n;
else
    [m,n]=size(A);
    A2  = A'*A;
end

%Compute Lipschitz constant
L = 2.1;   %for orthonormal A

%Convergence setup
MAX_ITER = 2000;
ABSTOL = 1e-7;

%Initialize
if nargin<6
	x = zeros(n,1);
else
    x = x0;
end
t = 1;
u = x;

out.et = [];out.e = [];
tic;

for iter = 1 : MAX_ITER

    xm1 = x;	

    if(isobject(A))
        v = x - (1/L)*(A'*(A*x) - Aty);
    else
        v = x - (1/L)*(A2*x - Aty);
    end

    x = shrinkage_Lq(v, q, lamda, L);  

   
    out.e  = [out.e norm(x-xtrue)/norm(xtrue)];
    out.et = [out.et toc];
        
    %Check for convergence
    if norm(x-xm1)<ABSTOL*sqrt(n)
        break;
    end

end


