function Z = nullsparsemacaulay(Z,problem,d,varargin)
    %NULLSPARSEMACAULAY   Sparse null space update of a Macaulay matrix.
    %   Z = NULLSPARSEMACAULAY(Z,problem,d) computes a numerical basis of 
    %   the null space of the Macaulay matrix of degree d that comprises 
    %   system, without building the corresponding Macaulay matrix.
    %
    %   Z = NULLSPARSEMACAULAY(...,tol) uses a user-specified tolerance for
    %   the rank decision.
    %
    %   Z = NULLSPARSEMACAULAY(...,basis) uses a specific monomial basis.
    %
    %   Z = NULLSPARSEMACAULAY(...,order) uses a specific monomial order.
    %
    %   See also NULL, MACAULAYUPDATE, NULLRECRMACAULAY, NULLITERMACAULAY.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Process the optional parameters:
    tol = 1e-10;
    basis = @monomial;
    order = @grevlex;
    for i = 1:nargin-3
        switch class(varargin{i})
            case "function_handle"
                if nargin(varargin{i}) == 1
                    order = varargin{i};
                else
                    basis = varargin{i};
                end
            case "double"
                tol = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end
    
    % Fix the important variables:
    [s,~,n,di,~,nRows,nCols] = properties(problem);

    % Compute the basis dilation:
    C = basis(ones(1,n),ones(1,n));
    bd = size(C,1);

    % Determine the number of shifts:
    shifts = zeros(s,1); % total number of shifts per matrix polynomial.
    newshifts = zeros(s,1); % number of new shifts per matrix polynomial.
    for k = 1:s
        shifts(k) = nbmonomials(d-di{k},n);
        newshifts(k) = shifts(k) - nbmonomials(d - di{k} - 1,n);
    end
    cumshifts = [0; cumsum(newshifts)];
    [q,nullity] = size(Z);
    
    % Create the monomial vectors:
    K = monomials(d,n);

    % Increase the basis matrix:    
    X = zeros(nRows*sum(newshifts),nullity);
    Y = zeros(nRows*sum(newshifts),nbmonomials(d,n)*nCols-q);
    for k = 1:s
        support = problem.supp{k};
        data = reshape(permute(problem.coef{k},[2 3 1]), ...
            [nRows size(support,1)*nCols]);
        S = basis(K(shifts(k)-newshifts(k)+1:shifts(k),:),support);
        c = reshape(kron(S(:,1).',ones(1,nCols)), ...
            [bd*size(support,1)*nCols, newshifts(k)]).';
        blocks = reshape(order(S(:,2:end)),[bd*size(support,1) ...
            newshifts(k)]).';
        col = kron((blocks-1)*nCols, ones(1,nCols)) + ...
            kron(ones(size(blocks)),1:nCols);
        colx = (col < q + 1);
        coly = (col > q);
        for l = 1:newshifts(k)
            cdata = c(l,:).*kron(ones(1,bd),data);
            shiftidx = cumshifts(k)+(l-1)*nRows+1:cumshifts(k)+l*nRows;
            X(shiftidx,:) = cdata(:,colx(l,:))*Z(col(l,colx(l,:)),:);
            Y(shiftidx,col(l,coly(l,:))-q) = cdata(:,coly(l,:));
        end
    end

    % Update the numerical basis matrix:
    [~,s,V] = svd([X Y]);
    r = nnz(s > tol);
    V = V(:,r+1:end);
    Z = [Z*V(1:nullity,:); V(nullity+1:end,:)];
end
