function pos = grlex(A)
    %GRLEX   Graded lexicographic monomial order.
    %   pos = GRLEX(A) returns the position of a monomial in the monomial 
    %   basis using the graded lexicographic order. The result is a k-by-1
    %   vector that contains the position for every n-variate monomial in
    %   the k-by-n matrix A. 
    %
    %   See also POSITION.
            
    % Copyright (c) 2024 - Christof Vermeersch
    % 
    % Note that the code is not yet vectorized.

    [nr,nc] = size(A);
    pos = zeros(nr,1);
    for k = 1:nr
        v = A(k,:);
        d = sum(v);
        order = nbmonomials(d,nc);
        for l = 1:nc
            order = order - nbmonomials(sum(v(l+1:nc))-1,nc-l);
        end
        pos(k) = order;
    end
end