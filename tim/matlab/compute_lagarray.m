function L = compute_lagarray(lagSet)

lagArrays = numel(lagSet);

% Check for errors and find out the
% maximum number of lags L.

lags = 1;
for i = 1 : lagArrays
    if ~isnumeric(lagSet{i})
        error('The members of LAGSET must be numeric arrays.');
    end
    iLags = numel(lagSet{i});
    if iLags > lags
        if lags > 1
            error(['The members of LAGSET must be either scalars ', ...
                'or arrays with the same number of elements.']);
        end
        lags = iLags;
    end
end

% Expand scalars to arrays of length L.

newLagSet = lagSet;

for i = 1 : lagArrays
    if numel(lagSet{i}) == 1
        newLagSet{i} = ones(1, lags) * lagSet{i};
    else
        newLagSet{i} = lagSet{i}(:)';
    end
end

% Form a numeric array from the cell-array.

L = cell2mat(newLagSet');
