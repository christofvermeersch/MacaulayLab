function [dgap,ma,gapsize] = gap(Z,n,varargin)
    %GAP   Gap zone of the basis matrix.
    %   dgap = GAP(Z,n) returns the first degree of the gap zone in the
    %   given basis matrix or NaN if there is no gap zone.
    %
    %   [dgap,ma,gapsize] = GAP(Z,n) also returns the number of affine
    %   solutions and the size of the gap zone (in number of degree
    %   blocks).
    %
    %   [...] = GAP(...,blocksize) considers the number of rows in the
    %   block rows (cf., the null space of block Macaulay matrix).
    %
    %   [...] = GAP(...,tol) uses a user-specified tolerance for the
    %   required rank checks.
    %
    %   [...] = GAP(...,isRecursive) selects the recursive or iterative 
    %   approach based on the flag. The default setting (flag = false) is 
    %   not to use the recursive approach.
    %
    %   See also STDMONOMIALS, COLUMNCOMPR, GAPSTDMONOMIALS.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    tol = 1e-10;
    isRecursive = false;
    blocksize = 1;
    for i = 1:nargin-2
        switch class(varargin{i})
            case 'double'
                if varargin{i} < 1
                    tol = varargin{i};
                else
                    blocksize = varargin{i};
                end
            case "logical"
                isRecursive = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end
    
    % Determine gap, number of affine solutions, and size of the gap:
    gapsize = 0;
    dgap = NaN;
    if isRecursive % recursive approach.
        q = size(Z,2);
        nm = blocksize*nbmonomials(0,n);
        Y = null(Z(1:nm,:),tol);
        ma = q - size(Y,2); 
        k = 1;
        while blocksize*nbmonomials(k,n) <= size(Z,1)
            nm1 = blocksize*nbmonomials(k-1,n)+1;
            nm2 = blocksize*nbmonomials(k,n);
            Y = nullrecrrow(Y,Z(nm1:nm2,:));
            r = q - size(Y,2);
            if ma == r
                gapsize = gapsize + 1;
                dgap = k - gapsize + 1;
            else
                if ~isnan(dgap)
                    break
                else
                    ma = r;
                end
            end
            k = k + 1;
        end
    else % iterative approach.
        ma = 0;
        k = 1;
        while blocksize*nbmonomials(k,n) <= size(Z,1)
            nm = blocksize*nbmonomials(k,n);
            r = rank(Z(1:nm,:),tol);
            if ma == r
                gapsize = gapsize + 1;
                dgap = k - gapsize + 1;
            else
                if ~isnan(dgap)
                    break
                else
                    ma = r;
                end
            end
            k = k + 1;
        end
    end
end
