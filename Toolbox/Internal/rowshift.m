function S = rowshift(d,n,idx,poly,varargin)
    %ROWSHIFT   Matrix of shifts.
    %   S = ROWSHIFT(d,n,idx,poly) creates a sparse shift matrix in which
    %   rows with given indices are shifted by a certain shift polynomial.
    %
    %   S = ROWSHIFT(...,l) considers problems with l-by-1 eigenvectors.
    %
    %   S = ROWSHIFT(...,basis) performs the shift in a user-specified
    %   monomial basis.
    %
    %   S = ROWSHIFT(...,order) performs the shift in a user-specified
    %   monomial order.
    %
    %   See also ROWRECOMB.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional arguments:
    l = 1;
    basis = @monomial;
    order = @grevlex;
    for i = 1:nargin-4
        switch class(varargin{i})
            case "function_handle"
                if nargin(varargin{i}) < 2
                    order = varargin{i};
                else
                    basis = varargin{i};
                end
            case "double"
                l = varargin{i};
            otherwise 
                error("Argument type is not recognized.")
        end
    end

    % Determine basis dilation:
    bd = size(basis(ones(1,n),ones(1,n)),1);

    % Create extended monomials matrix:
    K = [monomials(d,n,order,l) kron(ones(1,nbmonomials(d,n)),1:l)'];

    % Determine shifts of all the rows with given indices:
    shifts = basis(K(idx,1:end-1),poly(:,2:end));
    positions = order(shifts(:,2:end));

    % Construct the matrix of shifts:
    rowIdx = kron(1:length(idx),ones(1,size(poly,1)*bd)).';
    colIdx = (positions-1)*l+kron(K(idx,end),ones(size(poly,1)*bd,1));
    values = shifts(:,1).*kron(ones(bd*length(idx),1),poly(:,1));
    S = sparse(rowIdx,colIdx,values,length(idx),nbmonomials(d,n,l));
end
