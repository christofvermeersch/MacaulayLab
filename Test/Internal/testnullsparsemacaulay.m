classdef testnullsparsemacaulay < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function polynomialtest(testCase)
            system = randsystem(3,3,3);
            M = macaulay(system,10);
            Z = nullsparsemacaulay(null(full(M)),system,11);
            M = macaulay(system,11);
            testCase.verifyLessThanOrEqual(norm(M*Z),testCase.tolerance);
        end

        function meptest(testCase)
            mep = randmep(2,2,5,4);            
            M = macaulay(mep,10);
            
            Z = nullsparsemacaulay(null(full(M)),mep,11);
            M = macaulay(mep,11);
            testCase.verifyLessThanOrEqual(norm(M*Z),testCase.tolerance);
        end

        function basistest(testCase)
            system = randsystem(4,2,4);
            M = macaulay(system,10,@chebyshev);
            Z = nullsparsemacaulay(null(full(M)),system,11,@chebyshev);
            M = macaulay(system,11,@chebyshev);
            testCase.verifyLessThanOrEqual(norm(M*Z),testCase.tolerance);
        end

        function ordertest(testCase)
            system = randsystem(5,3,5);
            M = macaulay(system,10,@grlex);
            Z = nullsparsemacaulay(null(full(M)),system,11,@grlex);
            M = macaulay(system,11,@grlex);
            testCase.verifyLessThanOrEqual(norm(M*Z),testCase.tolerance);
        end
    end
end
