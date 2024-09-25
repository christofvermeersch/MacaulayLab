function [A,B] = multmapcolumn(M,K,G,idx,varargin)
    %MULTMAPCOLUMN   Generalized multiplication matrices (column space).
    %   [A,B] = MULTMAPCOLUMN(Z,K,G,idx) sets-up the generalized
    %   multiplication matrices for the shift polynomial in G. The matrix M
    %   contains (a selection of) the basis matrix of the column space
    %   and the matrix K contains the monomials that belong to each row of
    %   that M. All the columns mentioned in idx are shifted.
    %
    %   [A,B] = MULTMAPCOLUMN(Z,K,G,idx) can also work with multiple shift
    %   polynomials combined in a cell array. The resulting matrices are
    %   then also combined in cell arrays.
    %
    %   [A,B] = MULTMAPCOLUMN(...,blocksize) specifies the size of a 
    %   block row in M (otherwise it is supposed to be equal to 1). If
    %   block size is larger than 1, then the matrix K identifies the
    %   vector indices with every monomial.
    %
    %   [A,B] = MULTMAPCOLUMN(...,basis) uses a user-specified monomial 
    %   basis.
    %
    %   [A,B] = MULTMAPCOLUMN(...,order) supplements the function with a
    %   function handle for the order of the monomials in K.
    %
    %   See also MULTMAPNULL.
    
    % Copyright (c) 2024 - Christof Vermeersch

    warning('Code is not optimized for efficiency.')
    
    % Process the optional parameters:
    blocksize = 1;
    basis = @monomial;
    order = @grevlex;
    for i = 1:nargin-4
        switch class(varargin{i})
            case "function_handle"
                if nargin(varargin{i}) < 2
                    order = varargin{i}; % order has one input.
                else
                    basis = varargin{i}; % basis has two inputs.
                end
            case "double"
                blocksize = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end

    % Set the correct number of shift problems:
    if iscell(G)
        isCell = true;
    else
        isCell = false;
        G = {G};
    end
    m = length(G);
    n = size(G{1},2) - 1;

    % Compute the degree:
    d = max(sum(K,2));

    % Determine the selection of rows in the basis:
    if blocksize == 1
        selection = order(K);
    else
        selection = (order(K(:,1:end-1))-1)*blocksize + K(:,end);
    end
    
    % Define the generalized multiplication matrices:
    A = cell(m,1);
    B = cell(m,1);
    for i = 1:m
        % Create the shiftmatrix:
        Sg = rowshift(d,n,idx,G{i},blocksize,order,basis);
        B2 = [];
        for k = 1:length(idx)
            Sghere = Sg(k,:);
            Sghere(idx) = 0;
            if ~any(Sghere)
                B2 = [B2 k];
            end
        end
        Sg = Sg(B2,1:length(idx));
        
        % Reorder the Macaulay matrix:
        [Q,aff,hit] = rowrecomb(d,n,idx,G{i},blocksize,order,basis);
        P = Q(selection,selection);
        N = M/P;

        % Construct R33 and R34:
        [~,R] = qr(fliplr(N),0);
        R = fliplr(R);
        R33 = R(size(N,2)-aff-hit+1:size(N,2)-aff,aff+1:aff+hit);
        R34 = R(size(N,2)-aff-hit+1:size(N,2)-aff,1:aff);
        
        % Creat the generalized multiplication maps:
        A{i} = [Sg; -R34];
        B{i} = [eye(aff-hit) zeros(aff-hit,hit); zeros(hit,aff-hit) R33]; 
    end

    % Change output to array when input shift polynomial is an array:
    if isCell == false
        A = A{1};
        B = B{1};
    end
end
