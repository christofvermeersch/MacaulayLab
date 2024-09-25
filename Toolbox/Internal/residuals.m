function [maxres,res] = residuals(solutions,problem,basis)  
    %RESIDUALS   Solution residuals for a problem.
    %   maxres = RESIDUALS(solutions,problem) computes the maximum residual
    %   of the computed solutions for a given problem.
    %
    %   [...,res] = RESIDUALS(...) also returns the individual residual
    %   value of every solution.
    %
    %   [...] = RESIDUALS(...,basis) uses a user-specified monomial basis.
    %   The default basis is the standard monomial basis.
    %
    %   See also VALUE.

    % Copyright (c) 2024 - Christof Vermeersch

    % Set the default monomial basis:
    if nargin < 3
        basis = @monomial;
    end
    
    % Compute the residuals for all the solutions:
    nSolutions = size(solutions,1);
    res = zeros(nSolutions,1);
    for i = 1:nSolutions
        singularvalues = svd(value(solutions(i,:),problem,basis));
        res(i) = singularvalues(end);
    end
    
    % Determine the maximum residual:
    maxres = max(res);
end