classdef problemstruct
    %PROBLEMSTRUCT   Class for generic problems.
    %   obj = PROBLEMSTRUCT(coef,supp) creates a structure that contains
    %   all the information of the problem defined by the coefficients and
    %   support. Every matrix polynomial has one entry in the s-by-1 cells
    %   coef and supp. The coefficients or k-by-l coefficient matrices are
    %   given in a r-by-k-by-l tensor and every mode corresponds to one row
    %   in the r-by-n support matrix.
    %
    %   obj = PROBLEMSTRUCT(...,basis) also specifies the monomial basis
    %   of the problem.
    %
    %   PROBLEMSTRUCT properties:
    %       coef     - tensor with the coefficient matrices.
    %       supp     - support of the coefficient matrices.
    %       basis    - monomial basis.
    %       type     - type of problem (systemstruct or mepstruct).
    %       n        - number of variables/parameters.
    %       s        - number of (matrix) equations.
    %       k        - number of rows of the coefficient matrices.
    %       l        - number of columns of the coefficient matrices.
    %       di       - total degree of every 
    %       dmax     - maximum total degree in the support.
    %       nnze     - number of non-zero coefficients.
    %       ma       - number of affine solutions.
    %       mb       - total number of solutions.
    %       multi    - flag set to true when there are multiplicities.
    %       posdim   - flag for positive-dimensional solution sets.
    %       solved   - flag set to true when the problem is solved.
    %       appid    - related applications.
    %       bibid    - references in BibTeX.
    %       comments - additional comments about the problem.
    %
    %   PROBLEMSTRUCT methods:
    %       properties - returns the most important properties.
    %       disp       - displays a short description.
    %       info       - gives a more elaborate description.
    %       structure  - visualizes the problem structure (support).
    %
    %   See also SYSTEMSTRUCT, MEPSTRUCT.

    % Copyright (c) 2024 - Christof Vermeersch

    properties
        coef
        supp
        n
        s
        k
        l
        di
        dmax
        nnze
        basis
        ma
        mb
        multi
        posdim
        solved
        appid
        bibid
        comments
    end

    methods
        function obj = problemstruct(coef,supp,basis)
            %PROBLEMSTRUCT   Constructor for the problemstruct class.
            %   obj = PROBLEMSTRUCT(coef,supp) creates a structure that 
            %   contains all the information of the problem defined by the 
            %   coefficients and support.
            %
            %   obj = PROBLEMSTRUCT(...,basis) also specifies the monomial 
            %   basis of the problem.

            % Process the optional argument:
            if nargin < 3
                obj.basis = "unknown";
            else
                obj.basis = basis;
            end
            
            % Create a problemstruct:
            obj.coef = coef;
            obj.supp = supp;
            obj.n = size(supp{1},2);
            obj.s = numel(coef);
            obj.di = cell(obj.s,1);
            obj.nnze = cell(obj.s,1);
            for k = 1:obj.s
                obj.di{k} = max(sum(supp{k},2));
                obj.nnze{k} = nnz(coef{k});
            end   
            obj.dmax = max(horzcat(obj.di{:}));
            [~,obj.k,obj.l] = size(coef{1});
        end 

        function [s,dmax,n,di,nnze,k,l] = properties(obj)
            %PROPERTIES   Most important properties of the problem.
            %   [s,dmax,n,di,nnze,k,l] = PROPERTIES(problem) returns the 
            %   most important properties of the problem.

            s = obj.s;
            n = obj.n;
            di = obj.di;
            dmax = obj.dmax;
            nnze = obj.nnze;
            k = size(obj.coef{1},2);
            l = size(obj.coef{1},3);
        end

        function output = disp(obj)
            %DISP   Description of the problem.
            %   DISP(problem) displays a short description of the problem.
        
            switch class(obj)
                case "systemstruct"
                    if obj.s < 2
                        if obj.n < 2
                            text = sprintf("     univariate polynomial" ...
                                + " with degree %d \n",obj.dmax);
                        else
                            text = sprintf("     multivariate " ...
                                + "polynomial with total degree %d \n", ...
                                obj.dmax);
                        end
                    else 
                       if obj.n < 2
                            text = sprintf("     system of %d " ... 
                                + "univariate polynomials with maximum" ...
                                + " degree %d \n",obj.s,obj.dmax);
                        else
                            text = sprintf("     system of %d " ... 
                                + "multivariate polynomials with " ...
                                + "maximum total degree %d \n",obj.s, ...
                                obj.dmax);
                        end
                    end
                    
                case "mepstruct"
                    if obj.n < 2
                        text = sprintf("     %d-parameter eigenvalue " ...
                            + "problem with degree %d \n",obj.n, ...
                            obj.dmax);
                    else
                        text = sprintf("     %d-parameter eigenvalue " ...
                            + "problem with total degree %d \n", ...
                            obj.n, obj.dmax);
                    end
                otherwise
                    if obj.k > 1 || obj.l > 1
                        text = sprintf("     problem with %d matrix " ...
                            + "equations in %d variables \n",obj.s,obj.n);
                    else 
                        text = sprintf("     problem with %d equations" ...
                            + " in %d variables \n",obj.s,obj.n);
                    end
            end
            
            if nargout > 0
                output = text;
            else
                fprintf(text);
            end
        end

        function output = info(obj)
            %INFO   Information about the problem.
            %   INFO(problem) displays a more elaborate description of the
            %   problem that contains all the relevant information  

            text = disp(obj);
            text = text + structure(obj);
            text = text + sprintf("    sparsity factor is %.2f \n",sparsity(obj));

            if nargout > 0
                output = text;
            else
                fprintf(text);
            end
        end

        function output = structure(obj)  
            %STRUCTURE   Visualization of the problem structure.
            %   STRUCTURE(problem) visualizes the problem structure.

            if obj.k > 1 || obj.l > 1
                text = newline;
                for i = 1:obj.s
                    support = obj.supp{i};
                    text = text + sprintf("        (");
                    for j = 1:size(support,1)
                        text = text + sprintf("A");
                        for m = 1:obj.n
                            text = text + sprintf("%d",support(j,m));
                        end
                        text = text + sprintf(" ");
                        if sum(support(j,:)) > 0
                            for m = 1:obj.n
                                if support(j,m) > 0
                                    text = text + sprintf("\x03bb%d",m);
                                end
                                if support(j,m) > 1
                                    text = text + sprintf("^%d ", ...
                                        support(j,m));
                                end
                            end
                        end
                        if j < size(support,1)
                            text = text + sprintf("+ ");
                        end
                    end
                    text = text + sprintf(") z = 0 \n");
                end
                text = text + newline;
                text = text + sprintf("    coefficient matrices are " ...
                    + "%d x %d matrices \n",obj.k,obj.l);
            else 
                text = newline;
                for i = 1:obj.s
                    text = text + sprintf("        p%d(x1,...,x%d) = 0" ...
                        + "     <- d%d = %d \n",i,obj.n,i,obj.di{i});
                end
                text = text + newline;
            end

            if nargout > 0
                output = text;
            else
                fprintf(text);
            end
        end
    end
end