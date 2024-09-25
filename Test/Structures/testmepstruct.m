classdef testmepstruct < matlab.unittest.TestCase
    methods (Test)
        % Test the class constructor:
        function constructortest(testCase)
            A = randn(3,2);
            B = randn(3,2);
            C = randn(3,2);
            tensor = zeros(3,3,2);
            tensor(1,:,:) = A;
            tensor(2,:,:) = B;
            tensor(3,:,:) = C;
            mep = mepstruct({A,B,C},[0 0; 0 1; 1 0]);
            testCase.verifyClass(mep,'mepstruct');
            testCase.verifyInstanceOf(mep,'problemstruct');
            testCase.verifyEqual(mep.mat,{A,B,C});
            testCase.verifyEqual(mep.coef,{tensor});
            testCase.verifyEqual(mep.supp,{[0 0; 0 1; 1 0]});
            testCase.verifyEqual(mep.basis,"unknown");
        end
    end
end
