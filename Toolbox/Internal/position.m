function pos = position(A,order)
    %POSITION   Position of monomial.
    %   pos = POSITION(A) determines the position in the monomial basis of 
    %   the k monomials in the k-by-n matrix A. The function uses the 
    %   graded negative lexicographic order to assign a position. 
    %
    %   pos = POSITION(...,order) uses a user-specified monomial order.
    %   This function serves as an interface for the different implemented
    %   monomial orders.
    %
    %   See also GRNLEX, GREVLEX, GRLEX, GRVLEX.
            
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Set the default monomial order:
    if nargin < 2
        order = @grevlex;
    end
    
    % Determine the position:
    pos = order(A);
end
