function M = macaulay(problem,d,varargin)
    %MACAULAY   Construction of the (block) Macaulay matrix.
    %   M = MACAULAY(problem,d) constructs the (block) Macaulay matrix of
    %   a given degree in the standard monomial basis. The rows and columns
    %   are indexed by the monomial order from small to large.
    %
    %   M = MACAULAY(...,basis) uses a user-specified monomial basis.
    %
    %   M = MACAULAY(...,order) uses a user-specified monomial order.
    %
    %   See also MACAULAYUPDATE.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameter:
    basis = @monomial;
    order = @grevlex;
    for i = 1:nargin-2
        if strcmp(class(varargin{i}),"function_handle")
            if nargin(varargin{i}) == 1
                order = varargin{i};
            else
                basis = varargin{i};
            end
        else
            error("Argument type is not recognized.")
        end
    end
    
    % Fix the important variables:
    [s,~,n,di,~,nRows,nCols] = properties(problem);
    nElements = nRows*nCols;

    % Set the default basis (and compute its dilation):
    C = basis(ones(1,n),ones(1,n));
    bd = size(C,1);

    % Determine the number of shifts and elements:
    shifts = zeros(s,1); % number of shifts per equation.
    elements = zeros(s+1,1); % cumulative number of elements (including 0).
    for k = 1:s
        shifts(k) = nbmonomials(d - di{k},n);
        elements(k+1) = elements(k) ...
            + bd*nElements*size(problem.supp{k},1)*shifts(k);
    end
    cumshifts = [0; cumsum(shifts)];
    
    % Create the monomial vectors:
    K = monomials(d - min(horzcat(di{:})),n,order);
    
    % Compute the (block) Macaulay data:
    M = zeros(elements(end),3);
    rowIdx = 1:nRows;
    colIdx = 1:nCols;
    for k = 1:s
        support = problem.supp{k};
        data = reshape(problem.coef{k},[size(support,1) nElements]).';
        S = basis(K(1:shifts(k),:),support);
        blockRowIdx = ((cumshifts(k)+1:cumshifts(k+1)) - 1)*nRows;
        M1 = kron(blockRowIdx',ones(size(support,1)*bd*nElements,1)) ...
            + kron(ones(nCols*size(support,1)*bd*shifts(k),1),rowIdx');
        M2 = kron((order(S(:,2:end))-1)*nCols,ones(nElements,1)) ...
            + kron(ones(size(support,1)*bd*shifts(k),1), ...
            kron(colIdx',ones(nRows,1)));
        M3 = kron(S(:,1),ones(nElements,1)) ...
            .*kron(ones(shifts(k)*bd,1),data(:));
        M(elements(k)+1:elements(k+1),:) = [M1 M2 M3];
    end
    
    % Construct the (block) Macaulay matrix:
    nRows = sum(shifts)*nRows;
    nCols = nbmonomials(d,n)*nCols;
    M = sparse(M(:,1),M(:,2),M(:,3),nRows,nCols);
end
