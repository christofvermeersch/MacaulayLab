function m = bezout(system)
    %BEZOUT   Bézout number for a multivariate polynomial system.
    %   m = BEZOUT(system) computes the Bezout number for a 
    %   multivariate polynomial system.
    %
    %   See also KUSHNIRENKO, BKK.
    
    % Copyright (c) 2025 - Christof Vermeersch

    % Verify input requirements:
    if ~isa(system,"systemstruct") && ~isa(system,"problemstruct")
        error("Only systemstruct is allowed as input.");
    end
    
    % Compute Bézout number.
    m = prod(vertcat(system.di{:}));
end