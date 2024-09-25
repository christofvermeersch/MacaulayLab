classdef testoutputstruct < matlab.unittest.TestCase
    methods (Test)
       % Test the class constructor:
        function constructortest(testCase)
            time = ones(3,1);
            nullity = ones(20,1);
            shiftvalues = randn(3,1);
            details = outputstruct(time,nullity,shiftvalues);
            testCase.verifyClass(details,'outputstruct');
            testCase.verifyEqual(details.time,time);
            testCase.verifyEqual(details.nullity,nullity);
            testCase.verifyEqual(details.shiftvalues,shiftvalues);
        end
    end
end