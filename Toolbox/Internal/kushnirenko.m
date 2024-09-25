function m = kushnirenko(system)
    %KUSHNIRENKO   Solution bound for a multivariate polynomial system.
    %   m = KUSHNIRENKO(system) computes the Kushnirenko bound for a
    %   multivariate polynomial system.
    %       
    %   See also BEZOUT, BKK.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Verify input requirements:
    if ~isa(system,"systemstruct")
        error("Only systemstruct is allowed as input.");
    end

    % Determine the unique support of the system:
    uniqueSupport = unique(cat(1,system.supp{:}),"rows");

    % Compute the Kushnirenko solution bound:
    [~, volume] = convhulln(uniqueSupport);
    m = round(factorial(system.n)*volume);
end