classdef testgapstdmonomials < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function gapstdmonomials1test(testCase)
            c = [1 2 5 6];
            [dgap,ma] = gapstdmonomials(c,2,3);
            testCase.verifyEqual(ma,4)
            testCase.verifyEqual(dgap,NaN)
        end

        function gapstdmonomials2test(testCase)
            c = [1 2 5 6];
            [dgap,ma] = gapstdmonomials(c,3,3);
            testCase.verifyEqual(ma,4)
            testCase.verifyEqual(dgap,3)
        end

        function gapstdmonomials3test(testCase)
            c = [1 2 3 14];
            [dgap,ma] = gapstdmonomials(c,3,3);
            testCase.verifyEqual(ma,3)
            testCase.verifyEqual(dgap,2)
        end

        function gapstdmonomials4test(testCase)
            c = [1 2 3 30];
            [dgap,ma] = gapstdmonomials(c,3,3,2);
            testCase.verifyEqual(ma,3)
            testCase.verifyEqual(dgap,2)
        end
    end
end
