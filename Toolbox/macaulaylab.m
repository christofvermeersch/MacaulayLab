function [solutions,details] = macaulaylab(problem,dend,options)
    %MACAULAYLAB   Numerical solver that uses the Macaulay matrix.
    %   solutions = MACAULAYLAB(problem) solves the problem with default
    %   options. The input and output are given as a problemstruct and 
    %   solutionstruct, respectively.
    %
    %   solutions = MACAYLAYLAB(...,dend) stops the iterative procedure at
    %   the specified final degree.
    %   
    %   solutions = MACAULAYLAB(...,options) uses user-specified options.
    %   The options are given in key-value pairs. For example, solutions =
    %   MACAULAYLAB(...,verbose = true).
    %
    %   solutions = MACAULAYLAB(...,dend,options) combines both the final 
    %   degree and user-specified options, but since dend is a positional
    %   argument, it must be added before the option list.
    %
    %   [...,details] = MACAULAYLAB(...) also returns additional details
    %   about the algorithm and solutions.
    %
    %   MACAULAYLAB options:
    %       tol        - tolerance for all the rank decisions.
    %       clustertol - tolerance for the clustering step.
    %       polynomial - definition of the polynomial shift.
    %       basis      - function handle with monomial basis.
    %       order      - function handle with monomial order.
    %       subspace   - selected subspace.
    %       blocks     - approach to deal with the blocks.
    %       enlarge    - approach to enlarge the subspace.
    %       check      - approach to perform the rank checks. 
    %       clustering - flag for applying the clustering step.
    %       posdim     - flag for postive-dimensional solution sets.
    %       verbose    - flag for verbose output.
    %   
    %   MACAULAYLAB solutions fields:
    %       values - numerical values of the solutions.
    %       maxres - maximum residuals.
    %       res    - all residuals.
    %
    %   MACAULAYLAB details field:
    %       shiftvalues  - numerical values of the shift polynomial.
    %       nullity      - nullity in the different degrees.
    %       time         - time required to run the algorithm.
    %       multiplicity - flag set when multiplicities are detected.
    %       
    %   See also PROBLEMSTRUCT, SOLUTIONSTRUCT, OUTPUTSTRUCT.
    
    % Copyright (c) 2025 - Christof Vermeersch

    % Check the input and set the default options:
    arguments
        problem problemstruct
        dend double = 100
        options.tol double = 1e-10;
        options.clustertol double = 1e-3;
        options.polynomial double = [randn(problem.n,1) eye(problem.n)];
        options.basis function_handle = @monomial;
        options.order function_handle = @grevlex;
        options.subspace string = "null space";
        options.blocks string = "block-wise";
        options.enlarge string = "sparse";
        options.check string = "iterative";
        options.clustering logical = true;
        options.posdim logical = false;
        options.verbose logical = false;
    end

    % Parse text input:
    switch options.subspace
        case "null space"
            isNull = true;
        case "column space"
            isNull = false;
            options.enlarge = "recursive";
            warning("Approach to enlarge subspace is set to recursive.")
        otherwise
            error("Selected subspace invalid.")
    end
    switch options.blocks
        case "block-wise"
            isBlocks = true;
        case "row-wise"
            isBlocks = false;
        otherwise
            error("Blocks can only be block-wise or row-wise.")
    end
    switch options.check
        case "iterative"
            isRecursive = false;
        case "recursive"
            isRecursive = true;
        otherwise
            error("Check can only be iterative or recursive.")
    end

    % Fix the important variables:
    [~,dmax,n] = properties(problem);
    blocksize = size(problem.coef{1},3);

    % Display the properties of the problem (verbose only):
    if options.verbose
        print("about")
        print("problem",problem)
    end

    % Determine the degree of shift polynomial and start degree:
    shifts = cell(n+1,1);
    shifts{1} = options.polynomial;
    for i = 1:n
        shifts{i+1} = [1 zeros(1,i-1) 1 zeros(1,n-i)];
    end
    dpolynomial = max(sum(options.polynomial(:,2:end),2));
    dstart = max(dmax,dpolynomial);
    
    % Declare necessary variables and set stabilization flag:
    time = zeros(6,1);
    nullity = zeros(dend,1);
    gapsize = NaN;

    % Display the properties of the algorithm (verbose only):
    if options.verbose
        print("algorithm",options.enlarge,options.subspace,options.blocks);
    end

    % Build the initial subspace:
    tic
    M = full(macaulay(problem,dstart,options.basis,options.order));
    if isNull
        Z = null(M);
        nullity(dstart) = size(Z,2);
    else
        Z = NaN;
        r = rank(M);
        nullity(dstart) = size(M,2) - r;
    end
    time(1) = toc;

    % Display the properties of the algorithm (verbose only):
    if options.verbose
        print('iteration',dstart,nullity(dstart),0,gapsize)
    end

    for d = dstart+1:dend
        % Enlarge the subspace:
        tic
        [M,Z] = enlarge(M,Z,problem,d,isNull,options);
        time(1) = time(1) + toc;

        % Check the subspace:
        tic
        if isNull
            nullity(d) = size(Z,2);
        else
            r = rank(full(M));
            nullity(d) = size(M,2) - r;
        end
        time(2) = time(2) + toc;

        if options.posdim || nullity(d) == nullity(d-dpolynomial)
            % Look for a gap zone:
            tic 
            [dgap,ma,gapsize,idx,idxinf] = check(M,Z,d,n,blocksize, ...
                isNull,isBlocks,isRecursive,options);
            time(2) = time(2) + toc;

            % Break the iterative procedure when there is a gap zone:
            if gapsize >= dpolynomial
                if options.posdim
                    mb = Inf;
                else
                    mb = nullity(d);
                end
                if options.verbose
                    print('iteration',d,nullity(d),nullity(d-1),gapsize)
                end
                break
            end
        end

        % Display the properties of the algorithm (verbose only):
        if options.verbose
            print('iteration',d,nullity(d),nullity(d-1),gapsize)
        end

        % Throw error when no gap zone is found:
        if d == dend
            error('Gap zone is not found at the final degree.')
        end
    end

    % Perform a column compression (only for null space based approach):
    if isNull 
        tic
        Z = columncompr(Z,dgap + dpolynomial - 1,n,blocksize);
        time(3) = toc;
    end

    % Solve the problem:
    tic
    K = monomials(dgap + dpolynomial-1,n,blocksize);
    if blocksize > 1
        K = [K kron(ones(1,nbmonomials(dgap+dpolynomial-1,n)),1:blocksize)'];
    end
    switch options.subspace
        case "null space"
            if isnan(idx) % using a block-wise approach with all indices.
                idx = 1:blocksize*nbmonomials(dgap-1,n);
            end
            [A,B] = multmapnull(Z,K,shifts,idx,blocksize,options.basis, ...
                options.order);
        case "column space"
            [A,B] = multmapcolumn(Z,K,shifts,idx,idxinf,blocksize, ...
                options.basis,options.order);
        otherwise
            error("Solution subspace is not recognized.")
    end

    [Q,S] = schur(pinv(full(A{1}))*full(B{1}),"complex");
    shiftvalues = diag(S);
    values = zeros(ma,n);

    for i = 1:n
        S = Q'*pinv(A{i+1})*B{i+1}*Q;
        values(:,i) = diag(S);
    end
    time(4) = toc;

    % Perform the optional clustering.
    if options.clustering && ma > 1
        tic
        [clusteredvalues,labels] = multiplicity([shiftvalues values], ...
            options.clustertol);
        shiftvalues = clusteredvalues(:,1);
        values = clusteredvalues(:,2:end);
        if max(labels) < ma
            isMultiplicity = true;
        else
            isMultiplicity = false;
        end
        time(5) = toc;
    end
 
    tic
    [maxres,res] = residuals(values,problem,options.basis);
    time(6) = toc;

    % Store the solutions and output information:
    details = outputstruct(time,nullity,shiftvalues);
    if options.clustering && ma > 1
        details.multiplicity = isMultiplicity;
    end
    solutions = solutionstruct(values);
    solutions.maxres = maxres;
    solutions.res = res;
    solutions.ma = ma;
    solutions.mb = mb;

    % Display the properties of the solution candidates (verbose only):
    if options.verbose
        print("solutions",options.clustering,solutions,details)
    end
end

% Define all the required subfunctions:
function [M,Z] = enlarge(M,Z,problem,d,isNull,options)
    %ENLARGE   Next degree of the subspaces.
    %   [M,Z] = ENLARGE(M,Z,problem,d,isNull,options) results in a the next
    %   degree of the column and null space.

    if isNull % enlarge the null space.
        switch options.enlarge
            case 'iterative'
                M = full(macaulay(problem,d,options.basis,options.order));
                Z = null(M);
            case 'recursive'
                rows = size(M,1);
                M = macaulayupdate(M,problem,d,options.basis, ...
                    options.order);
                Z = nullrecrmacaulay(Z,full(M(rows+1:end,:)),options.tol);
            case "sparse"
                M = NaN;
                Z = nullsparsemacaulay(Z,problem,d,options.tol, ...
                    options.basis,options.order);
            otherwise
                error("Selected strategy to enlarge is not supported.")
        end
    else % enlarge the column space.
        switch options.enlarge
            case 'iterative'
                M = full(macaulay(problem,d,options.basis,options.order));
            case "recursive"
                M = macaulayupdate(M,problem,d,options.basis, ...
                    options.order);
            otherwise
                error("Selected strategy to enlarge is not supported.")
        end
        Z = NaN;
    end
end

function [dgap,ma,gapsize,idx,idxinf] = check(M,Z,d,n,vectorsize, ...
    isNull,isBlocks,isRecursive,options)
    %CHECK   Rank structure of the subspaces.
    %   [...] = CHECK(...) checks the rank structure of the subspaces.

    if isBlocks && isNull % find gap via blocked or non-blocked approach.
            % Set the standard monomials to NaN:
            idx = NaN;
            idxinf = NaN;
            [dgap,ma,gapsize] = gap(Z,n,vectorsize,isRecursive,...
                options.tol);
    else
        if isNull
            % Find the standard monomials:
            idx = stdmonomials(Z,options.tol);
        else
            % Find the standard monomials:
            sinv = stdmonomials(full(flip(M')),options.tol);
            idx = setdiff(1:size(M,2),size(M,2) - sinv + 1);
        end
        % Determine the number of affine solutions:
        [dgap,ma,gapsize] = gapstdmonomials(idx,d,n,vectorsize);
        idxinf = idx(ma+1:end);
        idx = idx(1:ma);
    end
end

function print(component,varargin)
    %PRINT   Output of MacaulayLab via fprintf.
    %   PRINT(component) generates the output of MacaulayLab for a certain 
    %   component of the solver. It is an auxiliary function that uses
    %   the built-in fprintf.
    %
    %   PRINT(...,varargin) also processes the additional information in
    %   the variable argument list.

    LINEWIDTH = 75;
    linestr = pad('', LINEWIDTH,'_');

    switch component
        case "about"
            clc
            fileId = fopen('about.txt');
            about = fscanf(fileId,'%c');
            fclose(fileId);
            % fprintf(linestr);
            fprintf("\n");
            fprintf(about);
            fprintf("\n");
            fprintf(linestr);
            fprintf("\n");

        case "problem"
            fprintf("\n");
            fprintf(pad("PROBLEM INFO",LINEWIDTH,"both"));
            fprintf("\n");
            info(varargin{1});
            fprintf(linestr);
            fprintf("\n"); 

        case "algorithm"
            fprintf("\n");
            fprintf(pad("ALGORITHM OUTPUT",LINEWIDTH,"both"));
            fprintf("\n");
            fprintf("    using a Macaulay matrix algorithm")
            fprintf(" with a %s construction of \n",varargin{1});
            fprintf("    the %s and %s shifts \n \n",varargin{2}, ...
                varargin{3});
            fprintf("    | degree        | nullity        | increase" ...
                + "       | gap zone      | \n")
            fprintf(['    |------------------------------------------' ...
            '-----------------------| \n'])
          
        case 'iteration'
            degree = sprintf("%d",varargin{1});
            degree = degree + pad('',13-strlength(degree),' ');
            nullity = sprintf("%d",varargin{2});
            nullity = nullity + pad('',14-strlength(nullity),' ');
            increase = sprintf("%d",varargin{2}-varargin{3});
            increase = increase + pad('',14-strlength(increase),' ');
            gapzone = sprintf("%d",varargin{4});
            gapzone = gapzone + pad('',13-strlength(gapzone),' ');
            fprintf('    | %s | %s | %s | %s | \n', ...
            degree, nullity, increase,gapzone)

        case "solutions"
            fprintf(linestr);
             fprintf("\n");
            fprintf("\n");
            fprintf(pad("SOLUTIONS",LINEWIDTH,"both"));
            fprintf("\n");
            info(varargin{2});
            if varargin{1}
                if varargin{3}.multiplicity
                    fprintf("    multiple solutions detected \n");
                else
                    fprintf("    no multiple solutions detected \n");
                end
            end
            fprintf("\n")
            fprintf("    total computation is %.4e seconds \n", ...
                sum(varargin{3}.time))
            fprintf("    maximum residual is %.4e \n",varargin{2}.maxres)
            fprintf(linestr);
            fprintf("\n"); 
    end 
end



