classdef testgap < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function gap1test(testCase)
            A = randn(1,1)*randn(1,5);
            B = randn(3,2)*randn(2,5);
            C = randn(10,3)*randn(3,5);
            Z = [A; B; A; B; A ; A; C];
            [dgap,ma] = gap(Z,3);
            testCase.verifyEqual(dgap,2)
            testCase.verifyEqual(ma,3);
        end

        function gap2test(testCase)
            A = randn(3,2)*randn(2,5);
            B = randn(9,2)*randn(2,5);
            C = randn(30,3)*randn(3,5);
            Z = [A; B; A; B; A ; A; C];
            [dgap,ma] = gap(Z,3,3,10e-10);
            testCase.verifyEqual(dgap,2)
            testCase.verifyEqual(ma,4);
        end

        function gap3test(testCase)
            A = randn(3,2)*randn(2,5);
            B = randn(9,2)*randn(2,5);
            C = randn(30,3)*randn(3,5);
            Z = [A; B; A + 10e-6*randn(3,5); B; A ; A; C];
            [dgap,ma] = gap(Z,3,3,10e-5);
            testCase.verifyEqual(dgap,2)
            testCase.verifyEqual(ma,4);
        end

        function gap4test(testCase)
            A = randn(3,2)*randn(2,5);
            B = randn(9,2)*randn(2,5);
            C = randn(30,3)*randn(3,5);
            Z = [A; B; A + 10e-6*randn(3,5); B; A ; A; C];
            [dgap,ma] = gap(Z,3,3,10e-10);
            testCase.verifyEqual(dgap,3)
            testCase.verifyEqual(ma,5);
        end
         
        function recursivetest(testCase)
            A = randn(3,2)*randn(2,5);
            B = randn(9,2)*randn(2,5);
            C = randn(30,3)*randn(3,5);
            Z = [A; B; A; B; A ; A; C];
            [dgap,ma] = gap(Z,3,3,true);
            testCase.verifyEqual(dgap,2)
            testCase.verifyEqual(ma,4);
        end

        function nantest(testCase)
            p1 = [1 1 0 0 0; 1 0 1 0 0; -1 0 0 0 0];
            p2 = [1 1 0 1 0; 1 0 1 0 1];
            p3 = [1 1 0 2 0; 1 0 1 0 2; -1 0 0 0 0];
            p4 = [1 1 0 3 0; 1 0 1 0 3];
            system = systemstruct({p1,p2,p3,p4});
            Z = null(full(macaulay(system,6)));
            dgap = gap(Z,4);
            testCase.verifyEqual(dgap,NaN);
        end

        function thesistest(testCase)
            p = [1 2 0; 1 1 1; -2 0 0];
            q = [1 0 2; 1 1 1; -2 0 0];
            system = systemstruct({p,q});
            Z = null(full(macaulay(system,4)));
            [dgap,ma] = gap(Z,2);
            testCase.verifyEqual(dgap,2);
            testCase.verifyEqual(ma,2);
        end
    end
end
