function s = sparsity(problem,zerotol)
    %SPARSITY   Sparsity factor of a problem.
    %   s = SPARSITY(problem) returns the sparsity factor of a problem,
    %   which is the number of non-zero elements divided by the possible
    %   number of elements if all coefficient (matrices) would be non-zero.
    %   
    %   s = SPARSITY(...,zerotol) ignores all coefficients smaller than the
    %   given tolerance.

    % Copyright (c) 2024 - Christof Vermeersch 

    % Process the optional argument:
    if nargin < 2
        zerotol = 1e-16;
    end

    % Compute the sparsity factor:
    num = 0;
    den = 0;
    for i = 1:problem.s
        nnze = nnz(abs(problem.coef{i}) >= zerotol);
        num = num + nnze;
        den = den + problem.k*problem.l* ...
            nbmonomials(problem.di{i},problem.n);
    end
    s = num/den;
end
