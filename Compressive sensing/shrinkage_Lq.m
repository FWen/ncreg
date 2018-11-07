function [x] = shrinkage_Lq(b, q, lam, L)

if (q==0)
    x = b;
    i1 = find(abs(b)<=sqrt(2*lam/L));
    x(i1) = 0;
elseif (q<1 && q>0)
    max_iter = 20;
    ABSTOL   = 1e-8;
    x    = zeros(length(b),1);
    beta = ( L / (lam*q*(1-q)) )^(1/(q-2));
    f1   = lam*q*beta.^(q-1) + L*beta - L*abs(b);
    i0   = find(f1<0);
    if ~isempty(i0) 
        b_u = abs(b(i0));
        x_u = b_u;
        for i=1:max_iter              
            deta_x = (lam*q*x_u.^(q-1) + L*x_u - L*b_u) ./ (lam*q*(q-1)*x_u.^(q-2) + L);
            x_u    = x_u - deta_x;
            if ( norm(deta_x) < sqrt(length(x_u))*ABSTOL )
                break;
            end
        end
        x_u = x_u .* sign(b(i0));
        i1 = find(L/2 * b_u.^2 - lam*abs(x_u).^q - L/2*(x_u-b(i0)).^2 < 0);
        x_u(i1) = 0;
        x(i0) = x_u;
    end 
elseif (q==1)
    x = sign(b) .* max(abs(b)-lam/L, 0);
end

end
