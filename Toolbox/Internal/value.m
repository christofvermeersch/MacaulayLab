function values = value(point,varargin)
    %VALUE   Evaluation of the (matrix) equation(s) in a certain point.
    %   values = VALUE(point,problem) computes the value(s) of the problem 
    %   in the given point. 
    %
    %   values = VALUE(point,coef,supp) also results in the value(s) of the
    %   problem in the given point, but the problem is given via its
    %   coefficients and support. The variables coef and supp are a tensor
    %   and a matrix, for a single (matrix) equation, or a cell array with
    %   tensors and matrices, for a system of (matrix) equations.
    %
    %   values = VALUE(...,basis) computes the value(s) of the 
    %   problem in the given point. The problem is given in a specified
    %   monomial basis, while otherwise the standard monomial basis is
    %   used.
    %
    %   See also RESIDUALS.

    % Copyright (c) 2024 - Christof Vermeersch
    %
    % Note that this function is not yet using the efficient implementation
    % of the evaluation algorithm.

    % Process the optional arguments:
    basis = @monomial;
    switch class(varargin{1})
        case {"problemstruct","systemstruct","mepstruct"}
            coef = varargin{1}.coef;
            supp = varargin{1}.supp;

            if nargin > 2
                basis = varargin{2};
            end
        case "cell"
            coef = varargin{1};
            supp = varargin{2};

            if nargin > 3
                basis = varargin{3};
            end
        case "double"
            coef = varargin(1);
            supp = varargin(2);

            if nargin > 3
                basis = varargin{3};
            end
        otherwise 
            error("Argument type is not recognized.")
    end

    % Translate the affine point to a row vector:
    [nr,nc] = size(point);
    if nr > nc
        point = point.';
    end

    % Evaluate all the (matrix) polynomials:
    s = length(coef);
    [~,k,l] = size(coef{1});
    values = zeros(s,k,l);
    for i = 1:s
        values(i,:,:) = matrixvalue(coef{i},supp{i},point,basis);
    end

    % Make this function compatible with both types of problems:
    values = squeeze(values);
end

% Define the required subfunctions:
function matrixvalues = matrixvalue(coef,supp,point,basis)
    %MATRIXVALUE   Evaluation of the (matrix) polynomial in a given point.
    %   matrixvalues = MATRIXVALUE(coef,supp,point,basis) computes the
    %   (matrix) values of the (matrix) polynomial defined by its
    %   coefficient tensor (coef) and support matrix (supp) in a point. The
    %   monomial basis is always supplemented to this function.

    [~,alpha,beta,phi0,phi1] = basis(point,point);
    pointbasis = basisvalue(max(sum(supp,2)), ...
        length(point),point,alpha,beta,phi0,phi1);
    
    matrixvalues = zeros(size(coef,2),size(coef,3));
    for k = 1:size(supp,1)
        matrixvalues = matrixvalues ...
            + squeeze(coef(k,:,:))*pointbasis(position(supp(k,:)));
    end
end

function pointbasis = basisvalue(d,n,point,alpha,beta,phi0,phi1)
    %BASISVALUE   Basis monomials evaluated in a point.
    %   pointbasis = BASISVALUE(d,n,point,alpha,beta,phi) returns the
    %   evaluation of every basis monomial in the point, based on the given
    %   basis recursion properties.

    univariate = zeros(d+1,n);
    for i = 1:n
        univariate(1,i) = point(i)^phi0(i+1);
        univariate(2,i) = point(i)^phi1(i+1);
        for j = 3:d+1
            univariate(j,i) = alpha(1)*point(i)^alpha(i+1) ...
                *univariate(j-1,i) ...
                + beta(1)*point(i)^beta(i+1)*univariate(j-2,i);
        end
    end
    
    K = monomials(d,n);
    pointbasis = ones(size(K,1),1);
    for i = 1:size(K,1)
        for j = 1:n
            pointbasis(i) = pointbasis(i) * univariate(K(i,j)+1,j);
        end
    end
end
