function m = shapiro(mep)
    %SHAPIRO   Solution bound for a multiparameter eigenvalue problem.
    %   m = SHAPIRO(mep) computes the Shapiro-Shapiro solution bound
    %   for a linear multiparameter eigenvalue problem.
    %
    %   See also HKP.

    % Copyright (c) 2024 - Christof Vermeersch

    % Verify input requirements:
    if ~isa(mep,"mepstruct")
        error("Only mepstruct is allowed as input.");
    end

    % Compute Shapiro-Shapiro number:
    m = nchoosek(mep.n+mep.l-1,mep.n);
end 
