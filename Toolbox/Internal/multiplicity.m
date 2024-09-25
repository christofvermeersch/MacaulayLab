function [X,E] = multiplicity(X,clustertol)
    %MULTIPLICITY   Clustered multiple solutions.
    %   X = MULTIPLICITY(X) uses a clustering algorithm to compare and
    %   combine multiple solutions. The solutions with multiplicity greater
    %   than one are replaced by their cluster centers.
    %
    %   [X,E] = MULTIPLICITY(X) also returns the cluster labels of every
    %   solution.
    %
    %   [...] = MULTIPLICITY(...,clustertol) uses a user-specified 
    %   clustering tolerance.
    %
    %   See also CLUSTER.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Process the optional argument:
    if nargin < 2
        clustertol = 1e-2;
    end
    
    % Cluster the solutions:
    Y = [real(X) imag(X)];
    E = cluster(linkage(Y),'cutoff', clustertol, 'criterion', 'distance');

    % Define the cluster centers for every cluster:
    nCluster = max(E);
    for i = 1:nCluster
        mask = (E==i);
        maskSize = nnz(mask);
        X(mask,:) = kron(mean(X(mask,:),1),ones(maskSize,1));
    end
end     