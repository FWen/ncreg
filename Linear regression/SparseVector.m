function [xs] = SparseVector(N,K)
    xs = zeros(N,1);
    kind = randperm(N);
    kind = kind(1:K);
    xs(kind) = randn(K,1);
    
    xs = xs/norm(xs);
end