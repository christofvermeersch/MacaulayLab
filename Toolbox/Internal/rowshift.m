function S = rowshift(d,n,rows,shift,varargin)
    %ROWSHIFT   Matrix of shifts.
    %   S = ROWSHIFT(d,n,rows,shift) creates a sparse shift matrix in which
    %   rows are shifted by a certain shift polynomial.
    %
    %   S = ROWSHIFT(...,blocksize) considers the block version instead.
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
    blocksize = 1;
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
                blocksize = varargin{i};
            otherwise 
                error("Argument type is not recognized.")
        end
    end

    % Create shift matrix:
    K = monomials(d,n,order);
    blockpos = 1:blocksize;
    Kext = [kron(K,ones(blocksize,1)) kron(ones(size(K,1),1),blockpos')];

    C = basis(ones(1,n),ones(1,n));
    bd = size(C,1);

    S = zeros(bd*length(rows)*size(shift,1),3);
    shifted = basis(Kext(rows,1:end-1),shift(:,2:end));
    S(:,1) = kron(1:length(rows),ones(1,size(shift,1)*bd)).';
    pos = order(shifted(:,2:end));
    S(:,2) = (pos-1)*blocksize+ kron(Kext(rows,end), ...
        ones(size(shift,1)*bd,1));
    S(:,3) = shifted(:,1).*kron(ones(bd*length(rows),1),shift(:,1));
    S = sparse(S(:,1),S(:,2),S(:,3),length(rows), ...
        nbmonomials(d,n,blocksize));
end
