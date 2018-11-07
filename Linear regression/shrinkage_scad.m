function [x] = shrinkage_scad(b, lam)

a = 3.7;
x = b;

i1 = find(abs(b)<=2*lam);
x(i1) = sign(b(i1)) .* max(abs(b(i1))-lam, 0);

i2 = find(2*lam<abs(b) & abs(b)<=a*lam);
x(i2) = ((a-1)*b(i2)-sign(b(i2))*a*lam)./(a-2);

end
