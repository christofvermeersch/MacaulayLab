function M = macaulayupdate(N,problem,d,varargin)
    %MACAULAYUPDATE   Update of the (block) Macaulay matrix.
    %   M = MACAULAYUPDATE(N,problem,d) constructs the (block) Macaulay
    %   of a higher degree, given the (block) Macaulay matrix at the
    %   previous degree.
    %
    %   M = MACAULAYUPDATE(...,basis) uses a user-specified monomial basis.
    %
    %   M = MACAULAYUPDATE(...,order) uses a user-specified monomial order.
    %
    %   See also MACAULAY.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    basis = @monomial;
    order = @grevlex;
    for i = 1:nargin-3
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
    
    % Determine the number of shifts:
    shifts = zeros(s,1); % Number of shifts per polynomial.
    newshifts = zeros(s,1);
    elements = zeros(s+1,1);
    for k = 1:s
        shifts(k) = nbmonomials(d-di{k},n);
        newshifts(k) = shifts(k) -  nbmonomials(d-di{k}-1, n);
        elements(k+1) = elements(k) ...
            + bd*nElements*size(problem.supp{k},1)*newshifts(k);
    end
    cumshifts = [0; cumsum(newshifts)];

    % Create the monomial vectors:
    K = monomials(d - min(horzcat(di{:})),n);
    
    % Compute the (block) Macaulay data:
    M = zeros(elements(end),3);
    rowIdx = 1:nRows;
    colIdx = 1:nCols;
    for k = 1:s
        support = problem.supp{k};
        data = reshape(problem.coef{k},[size(support,1) nElements]).';
        S = basis(K(shifts(k)-newshifts(k)+1:shifts(k),:),support);
        blockRowIdx = ((cumshifts(k)+1:cumshifts(k+1)) - 1)*nRows;
        M1 = kron(blockRowIdx',ones(size(support,1)*bd*nElements,1)) ...
            + kron(ones(nCols*size(support,1)*bd*newshifts(k),1),rowIdx');
        M2 = kron((order(S(:,2:end))-1)*nCols,ones(nElements,1)) ...
            + kron(ones(size(support,1)*bd*newshifts(k),1), ...
            kron(colIdx',ones(nRows,1)));
        M3 = kron(S(:,1),ones(nElements,1)) ...
            .*kron(ones(newshifts(k)*bd,1),data(:));
        M(elements(k)+1:elements(k+1),:) = [M1 M2 M3];
    end
    
    % Update the (block) Macaulay matrix:
    orows = size(N,1);
    ocols = size(N,2);
    nrows = sum(newshifts)*nRows;
    ncols = nbmonomials(d,n)*nCols;
    M = sparse(M(:,1),M(:,2),M(:,3),nrows,ncols);
    M = [N zeros(orows, ncols-ocols); M];
end
