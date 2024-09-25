function [A,B] = multmapnull(Z,K,G,idx,varargin)
    %MULTMAPNULL   Generalized multiplication matrices (right null space).
    %   [A,B] = MULTMAPULL(Z,K,G,idx) sets-up the generalized
    %   multiplication matrices for the shift polynomial in G. The matrix Z
    %   contains (a selection of) the basis matrix of the null space and
    %   the matrix K contains the monomials that belong to each row of that
    %   Z. All the rows mentioned in idx are shifted. 
    %
    %   [A,B] = MULTMAPULL(Z,K,G,idx) can also work with multiple shift
    %   polynomials combined in a cell array. The resulting matrices are
    %   then also combined in cell arrays.
    %
    %   [A,B] = MULTMAPULL(...,blocksize) specifies the size of a block
    %   row in Z (otherwise it is supposed to be equal to 1). If block size
    %   is larger than 1, then the matrix K identifies the vector indices
    %   with every monomial. 
    %
    %   [A,B] = MULTMAPULL(...,basis) uses a user-specified monomial
    %   basis.
    %
    %   [A,B] = MULTMAPULL(...,order) supplements the function with
    %   a function handle for the order of the monomials in K.
    %
    %   See also MULTMAPCOLUMN.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
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
        % Create first generalized multiplication map:
        A{i} = Z(idx,:);

        % Create second generalized multiplication map:
        S = rowshift(d,n,idx,G{i},blocksize,basis,order);
        R = S(:,selection);
        Y = zeros(length(idx),size(Z,2));
        for j = 1:length(idx)
            Y(j,:) = R(j,:)*Z;
        end
        B{i} = Y;
    end
   
    % Change output to array when input shift polynomial is an array:
    if isCell == false
        A = A{1};
        B = B{1};
    end
end
