function m = hkp(mep)
    %HKP   Solution bound for a multiparameter eigenvalue problem.
    %   m = HKP(mep) computes the Hochstenbach-Kosir-Plestenjak (HKP) 
    %   solution bound for a polynomial multiparameter eigenvalue problem.
    %
    %   See also SHAPIRO.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Verify input requirements:
    if ~isa(mep,"mepstruct")
        error("Only mepstruct is allowed as input.");
    end
    
    % Compute the HKP number:
    m = mep.dmax^mep.n*nchoosek(mep.n+mep.l-1,mep.n);
end
