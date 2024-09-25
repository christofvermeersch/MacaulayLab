function [C,alpha,beta,phi0,phi1] = legendre(A,B)
    %LEGENDRE   Multivariate product Legendre basis.
    %   C = LEGENDRE(A,B) computes the product of the two Legendre
    %   monomials in A and B (= 1-by-n vectors). If A and/or B are
    %   polynomials, then all monomials of B (= k-by-n matrix) are
    %   multiplied by every monomial of A (= l-by-n matrix) and stacked on
    %   top of each other in the polynomial C (= r-by-(n+1) matrix). The
    %   first column of C contains the corresponding coefficients of each
    %   computed monomial.
    %
    %   [...,alpha,beta,phi0,phi1] = LEGENDRE(...) also returns the  
    %   alpha(x) and beta(x) of the recursion relation for the Legendre 
    %   basis, i.e., P_{k + 1}(x) = alpha(x) * P_k(x) + beta(x) 
    %   P_{k - 1}(x). The first two basis functions are given in phi0(x)
    %   and phi1(x).
    %
    %   See also SHIFT.

    % Copyright (c) 2024 - Christof Vermeersch

    error('Basis operation is not yet implemented.')    
end    
