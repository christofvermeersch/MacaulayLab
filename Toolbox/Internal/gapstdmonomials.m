function [dgap,ma,gapsize] = gapstdmonomials(c,d,n,blocksize)
    %GAPSTDMONOMIALS   Gap zone via the standard monomials.
    %   dgap = GAPSTDMONOMIALS(c,d,n) returns the first degree of the gap
    %   zone in the list of standard monomials or NaN if there is no gap 
    %   zone. The monomials up to total degree d are considered.
    %
    %   [dgap,ma,gapsize] = GAPSTDMONOMIALS(c,d,n) also returns the 
    %   number of affine solutions and the size of the gap zone (in number 
    %   of degree blocks).
    %
    %   [...] = GAPSTDMONOMIALS(...,blocksize) considers the number of 
    %   rows in the block rows (cf., the null space of block Macaulay 
    %   matrix).
    %
    %   See also STDMONOMIALS, COLUMNCOMPR, GAP.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameter:
    if nargin < 4
        blocksize = 1;
    end
    
    % Determine gap, number of affine solutions, and size of the gap:
    ma = 0;
    gapsize = 0;
    dgap = NaN;
    for k = 0:d
        nm = blocksize*nbmonomials(k,n);
        r = sum(c<=nm);
        if ma == r
            gapsize = gapsize + 1;
            dgap = k - gapsize + 1;
        else
            if ~isnan(dgap)
                break
            else
                ma = r;
            end
        end
    end     
end