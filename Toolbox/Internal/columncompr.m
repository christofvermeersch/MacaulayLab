function Z = columncompr(Z,dcompr,n,varargin)
    %COLUMNCOMPR   Column compression.
    %   Z = COLUMNCOMPR(Z,dcompr,n) compresses a basis matrix of the
    %   (right) null space of a (block) Macaulay matrix with blocksize 
    %   equal to 1.
    %
    %   Z = COLUMNCOMPR(...,blocksize) considers a block size equal to the
    %   variable blocksize.
    %
    %   Z = COLUMNCOMPR(...,tol) uses a user-defined tolerance.

    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    blocksize = 1;
    tol = 1e-10;
    for i = 1:nargin-3
        if strcmp(class(varargin{i}),"double")  
            if varargin{i} < 1
                tol = varargin{i};
            else
                blocksize = varargin{i};
            end
        else
            error("Argument type is not recognized.")
        end
    end
    
    % Perform the column compression:
    nbrows = blocksize*nbmonomials(dcompr,n);
    [U,S,~] = svd(Z(1:nbrows,:));
    ma = nnz(diag(S) > tol);
    Z = U(:,1:ma);
end
