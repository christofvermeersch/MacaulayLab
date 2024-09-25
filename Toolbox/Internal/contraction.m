function q = contraction(p,d,n,varargin)
    %CONTRACTION   Matrix representation of a polynomial.
    %   q = CONTRACTION(p,d,n) contracts a polynomial p (with maximum
    %   total degree d and n variables) given by its vector representation 
    %   (in the graded reverse lexicographic ordering) into its 
    %   corresponding coefficient matrix q.
    %
    %   q = CONTRACTION(...,tol) only considers coefficients larger
    %   than the user-specified tolerance.
    %
    %   q = CONTRACTION(...,order) uses a user-specified monomial order.
    %
    %   See also EXPANSION.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    tol = 1e-10;
    order = @grevlex;
    for i = 1:nargin-3
        switch class(varargin{i})
            case "function_handle"
                order = varargin{i};
            case "double"
                tol = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end
       
    % Contract the polynomial:
    p = fliplr(p);
    mask = find(abs(p) > tol);
    K = monomials(d,n,order);
    q = [p(mask).' K(mask,:)];    
end
