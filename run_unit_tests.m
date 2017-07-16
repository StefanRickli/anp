%% define test files to be run

tests = {};

tests{end+1} = testsuite('unit_tests/anp_d_contour_test.m');
tests{end+1} = testsuite('unit_tests/anp_usage_examples_test.m');
tests{end+1} = testsuite('unit_tests/anp_mimo_test.m');
tests{end+1} = testsuite('unit_tests/anp_random_test.m');

%% init TestRunner
runner = matlab.unittest.TestRunner.withTextOutput;

%% run tests
res = {};
for ii = 1:length(tests)
    res{end+1} = runner.run(tests{ii});
end

%% display results
disp('-----------------------------------------------------------------------------------')

for ii = 1:length(res)
    disp(res{ii});
end