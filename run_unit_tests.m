tests = {};

tests{end+1} = testsuite('unit_tests/anp_d_contour_test.m');
tests{end+1} = testsuite('unit_tests/anp_usage_examples_test.m');
tests{end+1} = testsuite('unit_tests/anp_mimo_test.m');
tests{end+1} = testsuite('unit_tests/anp_random_test.m');

runner = matlab.unittest.TestRunner.withTextOutput;

res = {};
for ii = 1:length(tests)
    res{end+1} = runner.run(tests{ii});
end

disp('-----------------------------------------------------------------------------------')

for ii = 1:length(res)
    disp(res{ii});
end