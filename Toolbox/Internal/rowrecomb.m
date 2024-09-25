function [P,aff,hit,rem] = rowrecomb(d,n,rows,shift,varargin)
    %ROWRECOMB   Row combination matrix.
    %   P = ROWRECOMB(d,n,rows,shift) creates a sparse recombination matrix 
    %   for a problem in which given rows are shifted by a certain shift
    %   polynomial.
    %
    %   [...,aff,hit,rem] = ROWRECOMB(...) also returns the size of the
    %   three individual groups of indices in the recombination matrix.
    %
    %   P = ROWRECOMB(...,blocksize) considers the block version instead.
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

    % Create empty recombination matrix:
    p = nbmonomials(d,n)*blocksize;
    P = zeros(p,p);

    % Create B matrix:
    P(1:length(rows),rows) = eye(length(rows));
    rep = rows;

    Sg = rowshift(d,n,rows,shift,blocksize,order,basis);
    B2 = [];
    for k = 1:length(rows)
        Sghere = Sg(k,:);
        Sghere(rows) = 0;
        if any(Sghere)
            B2 = [B2 k];
        end
    end

    % Create C matrix:
    P(length(rows)+1:length(rows)+length(B2),:) = Sg(B2,:);
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
    P(length(rows)+length(B2)+1:end,notrep) = eye(length(notrep));  
    P = sparse(P);
    
    % Output additional information:
    aff = length(rows);
    hit = length(B2);
    rem = p-aff-hit;
end

