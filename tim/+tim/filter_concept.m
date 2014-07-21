function messageSet = filter_concept(data)

messageSet = {};

if ~isnumeric(data) 
	messageSet{end + 1} = 'must be a numeric array';
end

if mod(numel(data), 2) == 0
	messageSet{end + 1} = 'must have an odd number of elements';
end	

if sum(data(:)) == 0
	messageSet{end + 1} = 'must not sum to zero';
end	
