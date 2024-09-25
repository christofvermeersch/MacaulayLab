classdef testhkp < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function oneparametertest(testCase)
            actual = hkp(randmep(2,1,5,5));
            expected = 10;
            testCase.verifyEqual(actual,expected);
        end

        function twoparametertest(testCase)
            actual = hkp(randmep(2,2,3,2));
            expected = 12;
            testCase.verifyEqual(actual,expected);
        end

        function thesistest(testCase)
            A00 = [1 2; 3 4; 3 4];
            A10 = [2 1; 0 1; 1 3];
            A11 = [3 4; 2 1; 0 1];
            A02 = [1 2; 4 2; 2 1];
            supp = [0 0; 1 0; 0 2; 1 1];
            mep = mepstruct({A00,A10,A02,A11,},supp);
            actual = hkp(mep);
            expected = 12;
            testCase.verifyEqual(actual,expected);
        end
    end
end
