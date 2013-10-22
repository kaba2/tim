function test_generate(testName)

testDirectory = fileparts(mfilename('fullpath'));

i = 0;
inputName = '';
while true
    inputName = [testDirectory, '/test_', testName, '_in_', num2str(i), '.mat'];
    if exist(inputName, 'file') ~= 2
        % Found an input name which does not exist yet.
        break
    end
    i = i + 1;
end

className = ['timtest.test_', testName];
classFileName = [testDirectory, '/test_', testName, '.m'];
if exist(classFileName, 'file') ~= 2
    error(['Test file ', className, ' does not exist.']);
end

test = eval([className, '()']);
input = test.generate();
save(inputName, 'input');
