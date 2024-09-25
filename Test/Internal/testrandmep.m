classdef testrandmep < matlab.unittest.TestCase
    methods (Test)
        % Test the methods:
        function oneparametertest(testCase)
            mep = randmep(3,1,10,10);
            [s,dmax,n,di,nnze,k,l] = properties(mep);
            expected = {1,3,1,{3},{400},10,10};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
        end

        function twoparametertest(testCase)
            mep = randmep(2,2,5,4);
            [s,dmax,n,di,nnze,k,l] = properties(mep);
            expected = {1,2,2,{2},{120},5,4};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
        end
        
        function complextest(testCase)
            realmep = randmep(2,2,5,4,false);
            complexmep = randmep(2,2,5,4,true);
            defaultmep = randmep(2,2,5,4);
            testCase.verifyTrue(isreal(realmep.coef{1}));
            testCase.verifyFalse(isreal(complexmep.coef{1}));
            testCase.verifyTrue(isreal(defaultmep.coef{1}));
        end

        function basistest(testCase)
            defaultmep = randmep(2,2,5,4);
            chebyshevmep = randmep(2,2,5,4,@chebyshev);
            testCase.verifyEqual(defaultmep.basis,"unknown");
            testCase.verifyEqual(chebyshevmep.basis,@chebyshev);
        end 

        function sparsetest(testCase)
            supp = [0 0; 0 2; 2 0];
            mep = randmep(2,2,10,9,supp);
            [s,dmax,n,di,nnze,k,l] = properties(mep);
            expected = {1,2,2,{2},{270},10,9};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
            mat = mep.coef{1};
            testCase.verifyTrue(nnz(squeeze(mat(1,:,:))) == 90);
            testCase.verifyTrue(nnz(squeeze(mat(2,:,:))) == 0);
            testCase.verifyTrue(nnz(squeeze(mat(3,:,:))) == 0);
            testCase.verifyTrue(nnz(squeeze(mat(4,:,:))) == 90);
            testCase.verifyTrue(nnz(squeeze(mat(5,:,:))) == 0);
            testCase.verifyTrue(nnz(squeeze(mat(6,:,:))) == 90);
        end
    end
end
