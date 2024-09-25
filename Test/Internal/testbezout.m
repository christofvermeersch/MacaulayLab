classdef testbezout < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function univariatetest(testCase)
            actual = bezout(randsystem(1,5,1));
            expected = 5;
            testCase.verifyEqual(actual,expected);
        end

        function bivariatetest(testCase)
            actual = bezout(randsystem(2,4,2));
            expected = 16;
            testCase.verifyEqual(actual,expected);
        end

        function trivariatetest(testCase)
            actual = bezout(randsystem(3,3,3));
            expected = 27;
            testCase.verifyEqual(actual,expected);
        end

        function mixedtest(testCase)
            p = [1 5 0 1];
            q = [1 1 1 1];
            r = [1 0 2 3];
            actual = bezout(systemstruct({p,q,r}));
            expected = 90;
            testCase.verifyEqual(actual,expected);
        end

        function oldtest(testCase)
            testCase.verifyTrue(bezout(systemstruct({[1 5]})) == 5);
            testCase.verifyTrue(bezout(systemstruct({[1 3 0],[1 1 2]})) ...
                == 9);
            testCase.verifyTrue(bezout(systemstruct({[1 1 1 1], ...
                [1 0 3 0], [1 0 1 2]})) == 27);
        end

        function thesistest(testCase)
            p = [1 2 0; 1 0 2; -6 1 0; 7 0 0];
            q = [1 1 0; -1 0 1; -3 0 0];
            system = systemstruct({p,q});
            actual = bezout(system);
            expected = 2;
            testCase.verifyEqual(actual,expected);
        end
    end
end
