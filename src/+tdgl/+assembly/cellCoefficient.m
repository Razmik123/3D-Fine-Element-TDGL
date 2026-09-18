function values = cellCoefficient(coefficient, nCells, name)
%CELLCOEFFICIENT Expand a scalar or cellwise coefficient to a column.

if nargin < 3
    name = 'coefficient';
end
if isscalar(coefficient)
    values = repmat(double(coefficient), nCells, 1);
else
    values = double(coefficient(:));
end
if numel(values) ~= nCells || any(~isfinite(values))
    error('tdgl:assembly:InvalidCellCoefficient', ...
        '%s must be finite and scalar or contain one value per cell.', name);
end
end
