function [C,alpha,beta,phi0,phi1] = monomial(A,B)
    %MONOMIAL   Standard multivariate monomial basis.
    %   C = MONOMIAL(A,B) computes the product of the two monomials
    %   in A and B (= 1-by-n vectors). If A and/or B are polynomials, 
    %   then all monomials of B (= k-by-n matrix) are multiplied by every 
    %   monomial of A (= l-by-n matrix) and stacked on top of each other in
    %   the polynomial C (= r-by-(n+1) matrix). The first column of C
    %   contains the corresponding coefficients of each computed monomial.
    %
    %   [...,alpha,beta,phi0,phi1] = MONOMIAL(...) also returns the 
    %   alpha(x) and beta(x) of the recursion relation for the standard 
    %   monomial basis, i.e., x^{k + 1} = alpha(x) * x^k + beta(x) 
    %   x^{k - 1}. The first two basis functions are phi0(x) and phi1(x).
    %
    %   See also SHIFT.

    % Copyright (c) 2024 - Christof Vermeersch

    % Perform the multiplication:
    C = [ones(size(A,1)*size(B,1),1) ...
        kron(A,ones(size(B,1),1)) + kron(ones(size(A,1),1),B)];

    % Set the properties of the monomial basis:
    n = size(A,2);
    alpha = ones(1,n+1);
    beta = [0 zeros(1,n)];
    phi0 = [1 zeros(1,n)];
    phi1 = ones(1,n+1);
end
