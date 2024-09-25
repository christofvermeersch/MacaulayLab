function C = companion(poly)
    %COMPANION   Construction of the Companion matrix.
    %   C = COMPANION(p) constructs the companion matrix of a univariate
    %   polyomial in its matrix representation.
    %
    %   C = COMPANION(P) constructs the companion matrix of a univariate
    %   polyomial in its cell representation.
    %
    %   C = COMPANION(system) constructs the companion matrix of a
    %   univariate polynomial in its system representation.
    %
    %   See also ROOTS.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Transform the univariate polynomial into its system representation:
    switch class(poly)
        case "cell"
            poly = systemstruct(poly);
        case "double"
            poly = systemstruct({poly});
        case "systemstruct"
            % do nothing.
        otherwise
            error("Argument type is not recognized.")
    end

    % Check if the polynomial is univariate:
    if poly.n ~= 1
        error('The polynomial is not univariate!')
    end
    
    % Construct the companion matrix:
    matrix = [poly.coef{1} poly.supp{1}];
    p = fliplr(expansion(matrix));
    C = zeros(poly.dmax,poly.dmax);
    C(2:poly.dmax,1:poly.dmax-1) = eye(poly.dmax-1);
    C(:,poly.dmax) = -p(1:end-1).'./p(end);
end
