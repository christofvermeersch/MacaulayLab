function system = randsystem(s,dmax,n,varargin)
    %RANDSYSTEM   Dense system with random coefficients.
    %   system = RANDSYSTEM(s,dmax,n) generates a dense system with n
    %   variables and s multivariate polynomial equations (maximum total
    %   degree is dmax).
    %
    %   system = RANDSYSTEM(...,complex) adds the option for the 
    %   coefficients to be complex (= boolean).
    %
    %   system = RANDSYSTEM(...,supp) results in a multivariate polynomial
    %   system with only non-zero coeffients for the elements in the given
    %   support.
    %
    %   See also SYSTEMSTRUCT.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    basis = "unknown";
    isComplex = false;
    supp = cell(s,1);
    supp(:) = {monomials(dmax,n)};
    for i = 1:nargin-3
        switch class(varargin{i})
            case "function_handle"
                basis = varargin{i};
            case "logical"
                isComplex = varargin{i};
            case "cell"
                supp = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end

    % Create the systemstruct:
    eqs = cell(s,1);
    if ~isComplex % for real coefficients.
        for k = 1:s
            eqs{k} = [randn(size(supp{k},1),1) supp{k}];
        end
    else % for complex coefficients.
        for k = 1:s
            eqs{k} = [randn(size(supp{k},1),1) ...
                + 1i*randn(size(supp{k},1),1) supp{k}];
        end
    end
    system = systemstruct(eqs,basis);
end