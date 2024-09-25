function pos = grvlex(V)
    %GRVLEX   Graded inverse lexicographic monomial order.
    %   pos = GRVLEX(A) returns the position of a monomial in the monomial 
    %   basis using the graded inverse lexicographic order. The result is
    %   a k-by-1 vector that contains the position for every n-variate
    %   monomial in the k-by-n matrix A. 
    %
    %   See also POSITION.
            
    % Copyright (c) 2024 - Christof Vermeersch
    % 
    % Note that the code is not yet vectorized.
    
    [nr,nc] = size(V);
    pos = zeros(nr,1);
    for k = 1:nr
        v = V(k,:);
        d = sum(v);
        order = nbmonomials(d,nc);
        for l = nc:-1:1
            order = order - nbmonomials(sum(v(1:nc-l))-1,nc-l);
        end
        pos(k) = order;
    end
end