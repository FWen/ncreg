function [x,out] = scad_ista(A,y,lamda,xtrue,x0);
% solves
%
%   minimize || Ax - y ||_2^2 + \lambda || x ||_q^q
%
% inputs
%	A,y,lambda: CS variables
%	0<=q<=1
%	xtrue: for debug, for calculation of errors
%   x0: intialization
% outputs
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
MAX_ITER = 1000;
ABSTOL = 1e-7;

%Initialize
if nargin<5
	x = zeros(n,1);
else
    x = x0;
end
% t = 1;
u = x;

out.e = [];
tic;

for iter = 1 : MAX_ITER

    xm1 = x;	

    if(isobject(A))
        v = u - (1/L)*(A'*(A*u) - Aty);
    else
        v = u - (1/L)*(A2*u - Aty);
    end

    x = shrinkage_scad(v, lamda/L);

    u = x;
%     tp1 = (1 + sqrt(1+4*t^2))/2;
%     u = x + (t-1)/(tp1)*(x-xm1);
%     t = tp1;
    
%     out.e  = [out.e, norm(x-xtrue)/norm(xtrue)];
        
    %Check for convergence
    if norm(x-xm1)<ABSTOL*sqrt(n)
        break;
    end

end
