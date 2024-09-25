function pos = grevlex(A)
    %GREVLEX   Graded reverse lexicographic monomial order.
    %   pos = GREVLEX(A) returns the position of a monomial in the monomial 
    %   basis using the graded reverse lexicographic order. The result
    %   is a k-by-1 vector that contains the position for every n-variate
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
        for l = nc:-1:1
            order = order + nbmonomials(sum(v(1:nc-l+1))-1,nc-l+1);
        end
        pos(k) = order + 1;
    end
end