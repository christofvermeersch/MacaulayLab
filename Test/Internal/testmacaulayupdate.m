classdef testmacaulayupdate < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function systemtest(testCase)
            sys = randsystem(3,3,3);
            M1 = sort(macaulay(sys,5),1);
            M2 = sort(macaulayupdate(macaulay(sys,4),sys,5),1);
            testCase.verifyLessThanOrEqual(M1-M2,testCase.tolerance);
        end

        function meptest(testCase)
            mep = randmep(3,3,4,3);
            M1 = sort(macaulay(mep,5),1);
            M2 = sort(macaulayupdate(macaulay(mep,4),mep,5),1);
            testCase.verifyLessThanOrEqual(M1-M2,testCase.tolerance);
        end

        function basistest(testCase)
            sys = randsystem(4,3,3);
            M1 = sort(macaulay(sys,5,@chebyshev),1);
            M2 = sort(macaulayupdate(macaulay(sys,4,@chebyshev),sys,5, ...
                @chebyshev),1);
            testCase.verifyLessThanOrEqual(abs(M1-M2),testCase.tolerance);
        end

        function ordertest(testCase)
            mep = randmep(3,4,3,2);
            M1 = sort(macaulay(mep,4,@grnlex),1);
            M2 = sort(macaulayupdate(macaulay(mep,3,@grnlex),mep,4, ...
                @grnlex),1);
            testCase.verifyLessThanOrEqual(abs(M1-M2),testCase.tolerance);
        end
    end
end
