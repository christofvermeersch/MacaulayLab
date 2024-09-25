function Z = nullitermacaulay(problem,d,varargin)
    %NULLITERMACAULAY   Iterative null space basis matrix computation.
    %   Z = NULLITERMACAULAY(problem,d) computes the basis matrix of the
    %   null space of a Macaulay matrix of degree d for a given problem.
    %
    %   Z = NULLITERMACAULAY(...,approach) uses a specific approach to
    %   construct this basis matrix: "direct", "iterative", "recursive",
    %   and "sparse" are currently supported.
    %
    %   Z = NULLITERMACAULAY(...,tol) uses a user-specified tolerance for
    %   the rank decisions.
    %
    %   Z = NULLITERMACAULAY(...,basis) uses a user-specified monomial
    %   basis.
    %
    %   Z = NULLITERMACAULAY(...,order) uses a user-specified monomial
    %   order.
    %
    %   See also NULL, MACAULAYUPDATE, NULLRECRMACAULAY, NULLITERMACAULAY.

    % Copyright (c) 2024 - Christof Vermeersch

    % Process the optional parameters:
    basis = @monomial;
    order = @grevlex;
    tol = 1e-10;
    approach = "sparse";
    for i = 1:nargin-2
        switch class(varargin{i})
            case "function_handle"
                if nargin(varargin{i}) == 1
                    order = varargin{i};
                else
                    basis = varargin{i};
                end
            case "double"
                tol = varargin{i};
            case {"string","char"}
                approach = varargin{i};
            otherwise
                error("Argument type is not recognized.")
        end
    end

    % Construct a basis matrix for the null space:
    switch approach
        case "direct"
            M = macaulay(problem,d,basis,order);
            Z = null(full(M));

        case "iterative"
            M = macaulay(problem,problem.dmax,basis,order);
            Z = null(full(M));
            for i = problem.dmax+1:d
                M = macaulayupdate(M,problem,i,basis,order);
                Z = null(full(M));
            end

        case "recursive"
            M = macaulay(problem,problem.dmax,basis,order);
            Z = null(full(M));
            for i = problem.dmax+1:d
                M = macaulayupdate(M,problem,i,basis,order);
                Z = nullrecrmacaulay(Z,M,tol);
            end

        case "sparse"
            Z = null(full(macaulay(problem,problem.dmax,basis,order)));
            for i = problem.dmax+1:d
                Z = nullsparsemacaulay(Z,problem,i,tol,basis,order);
            end

        otherwise
            error("Approach to construct basis matrix is not recognized.")
    end
end
