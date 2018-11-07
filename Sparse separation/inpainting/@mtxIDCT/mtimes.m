function res = mtimes(A,x)

if A.adjoint == 0 %A*x
    res = idct(x);
else %At*x
    res = dct(x);
end
