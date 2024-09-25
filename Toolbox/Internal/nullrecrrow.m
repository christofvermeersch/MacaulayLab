function Z = nullrecrrow(Y,N,tol)
    %NULLRECRROW   Recursive null space update of a (block) row matrix.
    %   Z = NULLRECRROW(Y,N) builds a numerical basis Z of the null space 
    %   of a (block) row  matrix of degree d, using a numerical basis Y of 
    %   the null space of a lower degree (block) row matrix and the 
    %   difference N between both (block) row matrices.
    %
    %   Z = NULLRECRROW(...,tol) uses a user-specified tolerance to catch
    %   near-zero matrices.
    %
    %   See also NULL.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Process the optional parameter:
    if nargin < 3
        tol = 1e-10;
    end

    % Compute the updated basis matrix:
    X = N*Y;
    [~,S,V] = svd(X);
    r = nnz(diag(S) > tol);
    if r == 0
        V = eye(size(V,1));
    end
    Z = Y*V(:,r+1:end);
end