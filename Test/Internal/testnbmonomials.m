classdef testnbmonomials < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function fulltest(testCase)
            testCase.verifyEqual(nbmonomials(-1,3),0);
            testCase.verifyEqual(nbmonomials(0,3),1);
            testCase.verifyEqual(nbmonomials(1,3),4);
            testCase.verifyEqual(nbmonomials(2,2),6);
            testCase.verifyEqual(nbmonomials(3,3),20);
            testCase.verifyEqual(nbmonomials(4,4),70);
        end

        function blockedtest(testCase)
            testCase.verifyEqual(nbmonomials(3,3,3),60);
        end

        function thesistest(testCase)
            testCase.verifyEqual(nbmonomials(40,5),1221759);
        end
    end
end
