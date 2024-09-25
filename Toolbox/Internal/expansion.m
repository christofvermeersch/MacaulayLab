function q = expansion(p,varargin)
    %EXPANSION   Row vector representation of a polynomial.
    %   q = EXPANSION(p) expands a polynomial p given by its 
    %   coefficient matrix p into the corresponding row vector q (in the 
    %   graded reverse lexicographic ordering).
    %
    %   q = EXPANSION(...,tol) only considers coefficients larger
    %   than the user-specified tolerance.
    %
    %   q = EXPANSION(...,order) uses a user-specified monomial order.
    %   See also CONTRACTION.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    tol = 1e-10;
    order = @grevlex;
    for i = 1:nargin-1
        switch class(varargin{i})
            case "function_handle"
                order = varargin{i};
            case "double"
                tol = varargin{i};
            otherwise
                error("Argument type is not recognized.")
       end
    end

    % Create an empty expanded polynomial of appropriate length:
    n = size(p,2) - 1;
    d = max(sum(p(:,2:end),2));    
    q = zeros(1,nbmonomials(d,n));
    
    % Add non-zero coefficients to the expanded polynomial:
    mask = find(abs(p(:,1)) >= tol);
    coef = p(mask,1);
    pos = order(p(mask,2:end));
    q(end-pos+1) = coef;
end
