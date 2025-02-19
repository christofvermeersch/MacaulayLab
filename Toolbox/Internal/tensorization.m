function tensor = tensorization(mat)
    %TENSORIZATION   Conversion from matrix to tensor representation.
    %   tensor = TENSORIZATION(coef) creates a tensor of coefficient
    %   matrices from the cell array with coefficient matrices.
    % 
    %   See also MEPSTRUCT, PROBLEMSTRUCT.

    % Copyright (c) 2025 - Christof Vermeersch

    s = numel(mat);
    [nrows,ncols] = size(mat{1});
    tensor = permute(reshape(horzcat(mat{:}),[nrows,ncols,s]),[3,1,2]);
end