classdef testrandsystem < matlab.unittest.TestCase
    methods (Test)
        % Test the methods:
        function univariatetest(testCase)
            system = randsystem(1,4,1);
            [s,dmax,n,di,nnze,k,l] = properties(system);
            expected = {1,4,1,{4},{5},1,1};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
        end

        function trivariatetest(testCase)
            system = randsystem(3,4,3);
            [s,dmax,n,di,nnze,k,l] = properties(system);
            expected = {3,4,3,{4;4;4},{35;35;35},1,1};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
        end

        function overdeterminedtest(testCase)
            system = randsystem(3,3,2);
            [s,dmax,n,di,nnze,k,l] = properties(system);
            expected = {3,3,2,{3;3;3},{10;10;10},1,1};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
        end

        function complextest(testCase)
            realsystem = randsystem(1,3,2,false);
            complexsystem = randsystem(1,3,2,true);
            defaultsystem = randsystem(1,3,2);
            testCase.verifyTrue(isreal(realsystem.coef{1}));
            testCase.verifyFalse(isreal(complexsystem.coef{1}));
            testCase.verifyTrue(isreal(defaultsystem.coef{1}));
        end
            
        function basistest(testCase)
            defaultsystem = randsystem(8,8,8);
            chebyshevsystem = randsystem(8,8,8,@chebyshev);
            testCase.verifyEqual(defaultsystem.basis,"unknown");
            testCase.verifyEqual(chebyshevsystem.basis,@chebyshev);
        end

        function sparsetest(testCase)
            supp = {[0 0; 0 1; 2 0];[0 0; 3 0]};
            system = randsystem(2,3,2,supp);
            [s,dmax,n,di,nnze,k,l] = properties(system);
            expected = {2,3,2,{2;3},{3;2},1,1};
            testCase.verifyEqual({s,dmax,n,di,nnze,k,l},expected);
            testCase.verifyEqual(system.supp,supp);
        end
    end
end
