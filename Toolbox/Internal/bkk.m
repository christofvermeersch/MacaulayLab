function m = bkk(system)
    %BKK   Solution bound for a multivariate polynomial system.
    %   m = BKK(system) computes the Bernstein–Khovanskii–Kushnirenko (BKK)
    %   solution bound for a multivariate polynomial system.
    %
    %   See also BEZOUT, KUSHNIRENKO.
    
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % Note: this is a naive implementation of [1, Chapter 7, Theorem 4.12].
    %
    % [1] Cox, D. A., Little, J. B., and O'Shea, D. (2005) Using Algebraic
    % Geometry, vol. 185, Springer, New York, NY, USA.

    % Verify input requirements:
    if ~isa(system,"systemstruct")
        error("Only systemstruct is allowed as input.");
    end

    % Compute the BKK solution bound via Theorem 4.12:
    m = 0;
    for k = 1:system.n
        combinations = nchoosek(1:system.n,k);
        volume = 0;
        for l = 1:size(combinations,1)
            combo = combinations(l,:);
            B = cell(length(combo));
            for i = 1:length(combo)
                B{i} = system.supp{combo(i)};
            end
            [~, subvolume] = convhulln(minkowski(B));
            volume = volume + subvolume;
        end
        m = m + (-1)^(system.n-k)*volume;   
    end
    m = round(m);
end

% Define the required subfunction:
function sum = minkowski(A)
    %MINKOWSKI   Minkowski sum of vertices.
    %   sum = MINKOWSKI(A) computes the Minkowski sum of the vertices in A.

    sum = A{1};
    for k = 2:length(A)
        temp = [];
        for l = 1:size(A{k},1)
            temp = [temp; sum + A{k}(l,:)];
        end
        sum = unique(temp,"rows");
    end
end
