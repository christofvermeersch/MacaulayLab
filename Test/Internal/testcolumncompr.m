classdef testcolumncompr < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function systemtest(testCase)
            p = [1 2 0; 1 1 1; -2 0 0];
            q = [1 0 2; 1 1 1; -2 0 0];
            system = systemstruct({p,q});
            Z = null(full(macaulay(system,4)));
            W = columncompr(Z,3,2);
            testCase.verifySize(W,[10 3]);
            testCase.verifyEqual(rank(W),3);
        end

        function tolerancetest(testCase)
            Z = randn(10,3)*randn(3,5)+1e-7*randn(10,1)*randn(1,5);
            W = columncompr(Z,3,2,1e-4);
            testCase.verifySize(W,[10 3]);
            testCase.verifyEqual(rank(W),3);
        end

        function blockedtest(testCase)
            A00 = [1 2; 3 4; 3 4];
            A10 = [2 1; 0 1; 1 3];
            A11 = [3 4; 2 1; 0 1];
            A02 = [1 2; 4 2; 2 1];
            supp = [0 0; 1 0; 0 2; 1 1];
            mep = mepstruct({A00,A10,A02,A11},supp);
            Z = null(full(macaulay(mep,12)));
            W = columncompr(Z,3,2,2);
            testCase.verifySize(W,[20 9]);
            testCase.verifyEqual(rank(W),9);
        end
    end
end
