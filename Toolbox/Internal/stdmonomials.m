function v = stdmonomials(Z,tolerance)
    %STDMONOMIALS Standard monomials.
    %   v = STDMONOMIALS(Z) finds the standard monomials in the numerical 
    %   basis matrix Z of the (right) null space via rank checks.
    %
    %   v = STDMONOMIALS(...,tol) uses a user-specified tolerance for the
    %   required rank checks.
    %
    %   See also GAPSTDMONOMIALS.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameter:
    if nargin < 2
        tolerance = 1e-10; 
    end
    
    % Identify the standard monomials:
    ma = size(Z,2);
    vlist = zeros(size(Z,1),1);
    oldrank = 0;
    for k = 1:size(Z,1)
        newrank = rank(Z(1:k,:),tolerance);
        if newrank > oldrank
            vlist(k) = 1;
        end
        oldrank = newrank;
        if newrank == ma
            break
        end
    end
    v = find(vlist);
end
