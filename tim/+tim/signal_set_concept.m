function messageSet = SignalSet_concept(data)

messageSet = {};

if ~iscell(data)
	messageSet{end + 1} = 'must be a cell-array';
    return
end

signals = numel(data);

if signals == 0
	messageSet{end + 1} = 'must be non-empty';
    return
end

dimension = size(data{1}, 1);
maxDimension = 32;

for i = 1 : signals
	if size(data{i}, 1) ~= dimension
		messageSet{end + 1} = ...
            'must have equal dimensions for repeated trials';
	end

	if size(data{i}, 1) > maxDimension
		messageSet{end + 1} = ...
            ['must have signal-dimension at most ', ...
            int2str(maxDimension)];
	end
end

