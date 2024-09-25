classdef testmultmapcolumn < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (TestMethodSetup)
        function disablewarnings(~)
            warning('off','all')
        end
    end

    methods (TestMethodTeardown)
        function enablewarnings(~)
            warning('on','all')
        end
    end

    methods (Test)
        % Test the method:
        function blockedtest(testCase)
            A00 = [2 6; 4 5; 0 1];
            A10 = [1 0; 0 1; 1 1];
            A01 = [4 2; 0 8; 1 1];
            mep = mepstruct({A00,A01,A10},[0 0; 0 1; 1 0]);
            M = macaulay(mep,6);
            K = [monomials(6,2,2) kron(ones(1,nbmonomials(6,2)),1:2)'];
            [A,B] = multmapcolumn(M,K,[4 0 3],[1 2 3],2);
            actual = sort(eig(full(A),full(B)));
            expected = sort([-10.3979; 0.0007; -0.2928]);
            testCase.verifyLessThanOrEqual(actual-expected,1e-3);
        end

        function ordertest(testCase)
            p = [1 2 0; 1 0 2; -6 1 0; 7 0 0];
            q = [1 1 0; -1 0 1; -3 0 0];
            system = systemstruct({p,q});
            M = full(macaulay(system,3,@grnlex));
            K = monomials(3,2,@grnlex);
            [A,B] = multmapcolumn(M,K,[2 1 0; 3 0 1],[1 2],@grnlex);
            actual = sort(eig(full(A),full(B)));
            expected = sort([1; 11]);
            testCase.verifyLessThanOrEqual(actual-expected, ...
                testCase.tolerance);
        end

        function basistest(testCase)
            p = [1 2 0; 1 0 2; -6 1 0; 7 0 0];
            q = [1 1 0; -1 0 1; -3 0 0];
            system = systemstruct({p,q});
            M = full(macaulay(system,3,@chebyshev));
            K = monomials(3,2);
            [A,B] = multmapcolumn(M,K,[1 1 0],[1 2],@chebyshev);
            actual = sort(eig(full(A),full(B)));
            expected = sort([2.2500 + 0.8292i; 2.2500 - 0.8292i]);
            testCase.verifyLessThanOrEqual(actual-expected, ...
                1e-3);
        end

        function complextest(testCase)
            V = vandermonde(2,2,[2 -1i; 4i 1]);
            K = monomials(2,2);
            [actualA,actualB] = multmapnull(V,K,[2i 0 1], [1 2]);
            expectedA = [1 1; -1i 1];
            expectedB = [2 2i; -2i 2i];
            testCase.verifyEqual(actualA,expectedA);
            testCase.verifyEqual(actualB,expectedB);
        end

        function inputtest(testCase)
            K = monomials(5,3);
            M = macaulay(randsystem(3,2,3),5);
            [A1,B1] = multmapcolumn(M,K,ones(2,4),[1 2 4]);
            [A2,B2] = multmapcolumn(M,K,{ones(2,4)},[1 2 4]);
            [A3,B3] = multmapcolumn(M,K,{ones(2,4), ones(3,4)},[1 2 4]);
            testCase.verifyClass(A1,"double")
            testCase.verifyClass(B1,"double")
            testCase.verifyClass(A2,"cell");
            testCase.verifyClass(B2,"cell");
            testCase.verifyEqual(length(A2),1);
            testCase.verifyEqual(length(B2),1);
            testCase.verifyEqual(length(A3),2);
            testCase.verifyEqual(length(B3),2);
            testCase.verifyClass(A3{1},"double");
            testCase.verifyClass(B3{1},"double");
            testCase.verifyClass(A3,"cell");
            testCase.verifyClass(B3,"cell");
            testCase.verifyEqual(A1,A2{1});
            testCase.verifyEqual(B1,B2{1});
        end

        function usecasetest(testCase)
            p = [-1 2 0; 2 1 1; 1 0 2; 5 1 0; -3 0 1; -4 0 0];
            q = [1 2 0; 2 1 1; 1 0 2; -1 0 0];
            system = systemstruct({p,q});
            M = full(macaulay(system,3));
            K = monomials(3,2);
            [A,B] = multmapcolumn(M,K,[ones(2,1) eye(2)],[1 2 3 5]);
            actual = sort(eig(full(A),full(B)));
            expected = [-1; -1; 1; 1];
            testCase.verifyLessThanOrEqual(actual-expected, ...
                testCase.tolerance);
        end

        function secondusecasetest(testCase)
            p = [1 2 0; 1 0 2; -6 1 0; 7 0 0];
            q = [1 1 0; -1 0 1; -3 0 0];
            system = systemstruct({p,q});
            M = macaulay(system,3);
            [A,B] = multmapcolumn(M,monomials(3,2),[2 1 0; 3 0 1],[1 2]);
            actual = sort(eig(full(A),full(B)));
            expected = sort([1; 11]);
            testCase.verifyLessThanOrEqual(actual-expected, ...
                testCase.tolerance);
        end

        function thirdusecasetest(testCase)
            p = [1 2 0; 1 1 1; -2 0 0];
            q = [1 0 2; 1 1 1; -2 0 0];
            system = systemstruct({p,q});
            M = macaulay(system,5);
            K = monomials(5,2);
            remainder = setdiff(1:size(K,1),[11 16]);
            M = M(:,remainder);
            K = K(remainder,:);
            [A,B] = multmapcolumn(M,K,[1 1 0; 1 0 1],[1 2]);
            actual = sort(eig(full(A),full(B)));
            expected = sort([2; -2]);
            testCase.verifyLessThanOrEqual(actual-expected, ...
                testCase.tolerance);
        end
    end
end
