function C = shift(A,B,basis)
    %SHIFT   Result of shifting two monomials.
    %   C = SHIFT(A,B) computes the shift of the monomials/polynomial in B 
    %   by the monomials in A in the standard monomial basis. 
    %
    %   C = SHIFT(...,basis) considers the shift in a user-specified 
    %   monomial basis. This function serves as an interface for the 
    %   different implemented basis operations.
    %
    %   See also MONOMIAL, CHEBYSHEV, LEGENDRE.

    % Copyright (c) 2024 - Christof Vermeersch

    % Set the default monomial basis:
    if nargin < 3
        basis = @monomial;
    end

    % Compute the shifted monomial/polynomial:
    C = basis(A,B); % the function is an interface to the basis operation.
end
