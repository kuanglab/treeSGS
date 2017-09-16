function [ W ] = FusedLasso_SLEP( Y, X, W, gamma, lambda )
%FUSEDLASSO Summary of this function goes here
%   Detailed explanation goes here

max_iter = 10;

K = size(X, 2);

for iter = 1 : max_iter
    W_old = W;
    
    for k = 1 : K
        tmp = X(:, k)' * X(:, k);
        if tmp ~= 0
            W(k, :) = flsa_SLEP(X(:, k)' * (Y - X(:, [1:k-1 k+1:end]) * W([1:k-1 k+1:end], :)) / tmp, lambda / tmp, gamma / tmp);
        else
            W(k, :) = 0;
        end
    end
    
    cvg = norm(W - W_old, 'fro')/norm(W_old, 'fro');
%     disp(['FL:' num2str(iter) ' ' num2str(cvg)]);

    if cvg < 1e-3
        break;
    end
end

if iter == max_iter
    disp(['FusedLasso_SLEP didn''t converge in ' int2str(iter) ' iterations.']);
end

end