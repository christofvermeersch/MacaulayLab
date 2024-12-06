function P = rowrecomb(d,n,idx,poly,varargin)
    %ROWRECOMB   Row combination matrix.
    %   P = ROWRECOMB(d,n,idx,shift) creates a sparse recombination matrix 
    %   for a problem in which rows given by the indices are shifted by a 
    %   certain shift polynomial.
    %
    %   P = ROWRECOMB(...,l) considers problems with l-by-1 eigenvectors.
    %
    %   P = ROWRECOMB(...,basis) performs the shift in a user-specified
    %   monomial basis.
    %
    %   P = ROWRECOMB(...,order) performs the shift in a user-specified
    %   monomial order.
    %
    %   See also ROWSHIFT.
    
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

    % Create empty recombination matrix:
    p = nbmonomials(d,n)*l;
    P = zeros(p,p);

    % Create B matrix:
    P(1:length(idx),idx) = eye(length(idx));
    rep = idx;

    Sg = rowshift(d,n,idx,poly,l,order,basis);
    B2 = [];
    for k = 1:length(idx)
        Sghere = Sg(k,:);
        Sghere(idx) = 0;
        if any(Sghere)
            B2 = [B2 k];
        end
    end

    % Create C matrix:
    P(length(idx)+1:length(idx)+length(B2),:) = Sg(B2,:);
    for k = 1:length(B2)
        rowsinsg = find(Sg(B2(k),:));
        for l = 1:length(rowsinsg)
            if ~ismember(rowsinsg(l),rep)
                rep = [rep rowsinsg(l)];
                break
            end
        end
    end
    notrep = setdiff(1:p,rep);

    % Create D matrix:
    P(length(idx)+length(B2)+1:end,notrep) = eye(length(notrep));  
    P = sparse(P);
    
    % Output additional information:
    aff = length(idx);
    hit = length(B2);
    rem = p-aff-hit;
end
