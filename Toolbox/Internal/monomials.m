function K = monomials(d,n,varargin)
    %MONOMIALS   Matrix with all the monomial vectors.
    %   K = MONOMIALS(d,n) creates a matrix with all the monomial vectors
    %   of the monomial basis using the graded reverse lexicographic
    %   monomial order.
    %
    %   K = MONOMIALS(...,blocksize) repeats every monomial an amount of
    %   time equal to the blocksize.
    %
    %   K = MONOMIALS(...,order) uses a user-specified monomial order.
    %
    %   See also NBMONOMIALS.
        
    % Copyright (c) 2024 - Christof Vermeersch
    
    if d == 0
        K = zeros(1,n);
    else
        I = fliplr(eye(n,n));
        K = zeros(nbmonomials(d,n),n);
        K(2:n+1,:) = I;
        blockLength = ones(n,1);
        endPoint = 1;
        for di = 2:d
            blockLength = cumsum(blockLength,"reverse");
            blockPos = [0; cumsum(blockLength)];
            endPoint = endPoint + blockLength(1);
            for ni = n:-1:1
                K(endPoint+blockPos(ni)+1:endPoint+blockPos(ni+1),:) = ...
                    K(endPoint-blockLength(ni)+1:endPoint,:) + I(ni,:);             
            end
        end
    end

    % Order the monomials according to a specific monomial ordering.
    if nargin > 2
        blocksize = 1;
        for i = 1:nargin-2
            switch class(varargin{i})
                case "function_handle"
                    order = varargin{i};
                    K = K(order(K),:);
                case "double"
                    blocksize = varargin{i};
                otherwise
                    error("Argument type is not recognized.")
            end
        end
        K = kron(K,ones(blocksize,1));
    end
end
