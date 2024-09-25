function V = vandermonde(d,n,x,varargin)
    %VANDERMONDE   Multivariate (block) Vandermonde matrix.
    %   V = VANDERMONDE(d,n,x) constructs the multivariate Vandermonde
    %   matrix in the standard monomial basis and graded reverse
    %   lexicographic order. The Vandermonde matrix is built from the
    %   m different n-variate points in x (= m-by-n matrix).
    %
    %   V = VANDERMONDE(...,z) considers the blocked version instead, where
    %   the l-by-1 vectors are given in z (= m-by-l matrix). In the result,
    %   the elements of the vectors are also restructured as specified by
    %   the monomial order.
    %
    %   V = VANDERMONDE(...,order) works with a user-specified monomial
    %   order.
    %
    %   See also MONOMIALSMATRIX.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    z = ones(1,size(x,1)); % default is already transposed.
    order = @grevlex;
    isDifferentOrder = false;
    for i = 1:nargin-3
        switch class(varargin{i})
            case "function_handle"
                order = varargin{i};
                isDifferentOrder = true;
            case "double"
                z = varargin{i}.'; % transpose to put correctly in V.
            otherwise 
                error("Argument type is not recognized.")
        end
    end

    % Restructure the vectors:
    idx = order(eye(size(z,1)))-1;
    z = z(idx,:);

    % Construct Vandermonde matrix:
    K = monomials(d,n); 
    V = zeros(size(K,1)*size(z,1),size(x,1));
    for k = 1:size(x,1)
        V(:,k) = kron(prod(x(k,:).^K,2),z(:,k));
    end
    
    % Change the order of the rows when a different monomial order is used:
    if isDifferentOrder
        idx = kron(order(K)-1,size(z,1)*ones(size(z,1),1)) ...
            + kron(ones(1,size(K,1)),1:size(z,1))';
        V = V(idx,:);
    end
end