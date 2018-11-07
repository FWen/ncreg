function [x1,x2,out] = lq_lq_admm(A1,A2,y,mu,q1,q2,x01,x02,xtrue);
% lq_lq_admm solves
%
%   minimize \mu*|| x1 ||_q1^q1 + || x2 ||_q2^q2    s.t.   A1*x1 + A2*x2 - y  = 0
%
% Inputs
%	A1,A2,y,mu: 
%	1=>q1,q2=>0
%	xtrue: for debug, for calculation of errors
%   x1,x2: intialization
% Outputs
%	x1,x2: the recovery
%	out.e: the error with respect to the true
%	out.et: time index

%Convergence setup
MAX_ITER = 5e2; 
ABSTOL = 1e-7;

[m, L] = size(y);

rho  = 0.2;

if(isobject(A1))
    n1 = A1.n;
else
    n1 = size(A1,2);
end

if(isobject(A2))
    n2 = A2.n;
else
    n2 = size(A2,2);
end

%Initialize
if nargin<7
	x1 = zeros(n1,L);
    x2 = zeros(n2,L);
else
    x1 = x01;
    x2 = x02;
end

isquiet = 0;
if nargin<9
    isquiet = 1;
end

w = zeros(n2,L);

c1 = 2.1*norm(A1'*A1);
c2 = 2.1*norm(A2'*A2);

out.et = [];out.e = [];
tic;

for iter = 1:MAX_ITER

    x1m1 = x1;
    x2m1 = x2;
	    
    % x1-update         
    b1 = x1 - 1/c1*(A1'*(A1*x1 + A2*x2 - y - w/rho));
    x1 = repmat(prox_Lq_m(b1,q1,c1*rho/mu),[1,L]) .* b1;
    
    % x2-update
    b2 = x2 - 1/c2*(A2'*(A1*x1 + A2*x2 - y - w/rho));
    x2 = repmat(prox_Lq_m(b2, q2, c2*rho),[1,L]) .* b2; 
   
 	% w-update
 	w = w - rho*(A1*x1 + A2*x2 - y);
    
    if ~isquiet
        out.e  = [out.e norm(x1-xtrue)/norm(xtrue)];
        out.et = [out.et, toc];
    end
        
    %Check for convergence: when both primal and dual residuals become small
    if (norm(A1*x1 + A2*x2 - y) < sqrt(n1)*ABSTOL) 
        break;           
    end    
    
end

end
