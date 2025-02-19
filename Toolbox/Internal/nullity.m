function n = nullity(M,tol)
    %NULLITY   Nullity of a matrix.
    %   n = NULLITY(M) returns the nullity of the matrix M. 
    %
    %   n = NULLITY(...,tol) uses a tolerance for the rank determination. 
    %
    %   See also RANK.

    % Copyright (c) 2025 - Christof Vermeersch

    % Determine the rank
    if nargin < 2
        r = rank(M);
    else
        r = rank(M,tol);
    end
    
    % Determine the nullity:
    n = size(M,2) - r;
end