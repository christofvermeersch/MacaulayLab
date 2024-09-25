classdef testnullrecrmacaulay < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function polynomialtest(testCase)
            system = randsystem(3,3,3);
            M1 = macaulay(system,10);
            M2 = macaulayupdate(M1,system,11);
            
            Z = nullrecrmacaulay(null(full(M1)), ...
                full(M2(size(M1,1)+1:end,:)));
            testCase.verifyLessThanOrEqual(norm(M2*Z),testCase.tolerance);
        end

        function meptest(testCase)
            mep = randmep(2,2,5,4);            
            M1 = macaulay(mep,10);
            M2 = macaulayupdate(M1,mep,11);
            
            Z = nullrecrmacaulay(null(full(M1)), ...
                full(M2(size(M1,1)+1:end,:)));
            testCase.verifyLessThanOrEqual(norm(M2*Z),testCase.tolerance);
        end
    end
end
