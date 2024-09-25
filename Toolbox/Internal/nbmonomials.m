function s = nbmonomials(d,n,blocksize)
    %NBMONOMIALS   Number of monomials.
    %   s = NBMONOMIALS(d,n) returns the number of monomials in the 
    %   monomial basis for n variables and maximum total degree d.
    %
    %   s = NBMONOMIALS(...,blocksize) takes into account the blocksize.
    %
    %   See also MONOMIALS.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Compute the number of monomials:
    if d < 0
        s = 0; % avoid errors with negative degrees
    else
        s = round(exp(gammaln(d+n+1)-gammaln(n+1)-gammaln(d+1)));
    end

    % Take into account the blocksize:
    if nargin > 2
        s = s*blocksize;
    end
end
