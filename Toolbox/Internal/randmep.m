function mep = randmep(dmax,n,k,l,varargin)
    %RANDMEP   Dense multiparameter eigenvalue problem.
    %   mep = RANDMEP(dmax,n,k,l) generates a dense multiparameter 
    %   eigenvalue problem with n eigenvalues and random k-by-l coefficient
    %   matrices (maximum total degree is equal to dmax).
    %
    %   mep = RANDMEP(...,complex) adds the option for the coefficient
    %   matrices to be complex (= boolean).
    %
    %   mep = RANDMEP(...,supp) results in a multiparameter eigenvalue
    %   problem with only non-zero coefficient matrices for the elements in
    %   the given support.
    %
    %   See also MEPSTRUCT.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    basis = "unknown";
    isComplex = false;
    supp = monomials(dmax,n);
    for i = 1:nargin-4
        switch class(varargin{i})
            case "function_handle"
                basis = varargin{i};
            case "logical"
                isComplex = varargin{i};
            case "double"
                supp = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end

    % Create the mepstruct:
    dmax = max(sum(supp,2));
    s = nbmonomials(dmax,n);
    mat = cell(s,1);
    mat(:) = {zeros(k,l)};
    if ~isComplex % for real coefficient matrices.
        for i = grevlex(supp)'
            mat{i} = randn(k,l);
        end
    else % for complex coefficient matrices.
        for i = grevlex(supp)'
            mat{i} = randn(k,l) + 1i*randn(k,l);
        end
    end
    mep = mepstruct(mat,supp,basis);
end