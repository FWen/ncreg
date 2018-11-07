function [x] = prox_Lq_m(b, q, phi)

if (q==0)
    x = ones(size(b,1),1);
    i1 = find(1<=sqrt(2/phi./sum(b.^2,2)));
    x(i1) = 0;
elseif (q<1 && q>0)
    max_iter = 10;
    ABSTOL   = 1e-5;
    
    L    = phi*sum(b.^2,2).^(1-q/2);
    x    = zeros(length(b),1);
    beta = ( 2*(1-q)./L ).^(1/(2-q));
    tao  = beta + q*beta.^(q-1)./L;
    i0   = find(tao<=1);
    n    = length(i0);
    if n>0 
        b_u = ones(n,1);
        x_u = b_u;
        L2 = L(i0);
        for k=1:max_iter              
            deta_x = (q*x_u.^(q-1) + L2.*(x_u - b_u)) ./ (q*(q-1)*x_u.^(q-2) + L2);
            x_u    = x_u - deta_x;
            if (k>2 && norm(deta_x) < sqrt(n)*ABSTOL )
                break;
            end
        end
        x(i0) = x_u;
    end 
elseif (q==1)
    x = max(1-1./(phi*sqrt(sum(b.^2,2))), 0);
end

end
