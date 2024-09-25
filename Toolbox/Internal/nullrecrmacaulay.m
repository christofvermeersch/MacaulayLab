function Z = nullrecrmacaulay(Y,N,tol)
    %NULLRECRMACAULAY   Recursive null space update of a Macaulay matrix.
    %   Z = NULLRECRMACAULAY(Y,N) builds a numerical basis of the 
    %   null space  of the Macaulay  matrix of degree d, using a 
    %   numerical basis of the null space Y of a lower degree Macaulay 
    %   matrix and the difference N between both Macaulay matrices.
    %
    %   Z = NULLRECRMACAULAY(...,tol) uses a user-specified tolerance for
    %   the rank decision.
    %
    %   See also NULL, MACAULAYUPDATE, NULLSPARSEMACAULAY, 
    %   NULLITERMACAULAY.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Process the optional parameter:
    if nargin < 3
        tol = 10e-10;
    end

    % Compute the updated basis matrix:
    [q,n] = size(Y); 
    X = N(:,1:q)*Y;
    [~,s,V] = svd(full([X N(:,q+1:end)]));
    r = nnz(s > tol);
    V = V(:,r+1:end);
    Z = [Y*V(1:n,:); V(n+1:end,:)];
end

