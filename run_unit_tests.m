tests_d_contour =       testsuite('unit_tests/anp_d_contour_test.m');
tests_usage_examples =  testsuite('unit_tests/anp_usage_examples_test.m');

runner = matlab.unittest.TestRunner.withTextOutput;

res1 = runner.run(tests_d_contour);
res2 = runner.run(tests_usage_examples);

disp('-----------------------------------------------------------------------------------')
disp(res1);
disp(res2);
