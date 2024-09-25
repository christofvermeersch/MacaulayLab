function pos = grnlex(A)
    %GRNLEX   Graded negative lexicographic monomial order.
    %   pos = GRNLEX(A) returns the position of a monomial in the monomial 
    %   basis using the graded negative lexicographic order. The result is
    %   a k-by-1 vector that contains the position for every n-variate
    %   monomial in the k-by-n matrix A. 
    %
    %   See also POSITION.
            
    % Copyright (c) 2024 - Christof Vermeersch
    % 
    % Note that the code is not yet vectorized.
    
    [nr,nc] = size(A);
    pos = zeros(nr,1);
    for k = 1:nr
        v = A(k,:);
        order = 0;
        for l = 1:nc-1
            order = order + nbmonomials(sum(v(l:nc))-1,nc-l+1);
        end
        pos(k) = order + v(nc) + 1;
    end
end