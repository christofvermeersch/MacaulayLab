classdef testnullitermacaulay < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function directtest(testCase)
            system = randsystem(2,2,2);
            Z = nullitermacaulay(system,5,"direct");
            M = full(macaulay(system,5));
            testCase.verifyLessThanOrEqual(M*Z,testCase.tolerance);
        end

        function iterativetest(testCase)
            system = randsystem(2,2,2);
            Z = nullitermacaulay(system,5,'iterative');
            M = full(macaulay(system,5));
            testCase.verifyLessThanOrEqual(M*Z,testCase.tolerance);
        end

        function recursivetest(testCase)
            system = randsystem(2,2,2);
            Z = nullitermacaulay(system,5,'recursive');
            M = full(macaulay(system,5));
            testCase.verifyLessThanOrEqual(M*Z,testCase.tolerance);
        end

        function sparsetest(testCase)
            system = randsystem(2,2,2);
            Z = nullitermacaulay(system,5,'sparse');
            M = full(macaulay(system,5));
            testCase.verifyLessThanOrEqual(M*Z,testCase.tolerance);
        end

        function optionstest(testCase)
            mep = randmep(2,2,4,3);
            Z = nullitermacaulay(mep,5,"direct",@grnlex,@chebyshev);
            M = macaulay(mep,5,@grnlex,@chebyshev);
            testCase.verifyLessThanOrEqual(norm(M*Z),testCase.tolerance);
        end
    end
end
