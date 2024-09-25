classdef testmultmapnull < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function simpletest(testCase)
            V = vandermonde(2,2,[2 -1; 4 1]);
            K = monomials(2,2);
            [actualA,actualB] = multmapnull(V,K,[2 0 1; 3 1 0], [1 2]);
            expectedA = [1 1; -1 1];
            expectedB = [4 14; -4 14];
            testCase.verifyEqual(actualA,expectedA);
            testCase.verifyEqual(actualB,expectedB);
        end
        
        function blockedtest(testCase)
            V = vandermonde(2,2,[2 3; 4 5; 0 0],[2 1; 1 1; 1 0]);
            K = [monomials(2,2,2) kron(ones(1,nbmonomials(2,2)),1:2)'];
            [actualA,actualB] = multmapnull(V,K,[1 0 1],[1 2 3],2);
            expectedA = [1 1 0; 2 1 1; 3 5 0];
            expectedB = [3 5 0; 6 5 0; 9 25 0];
            testCase.verifyEqual(actualA,expectedA);
            testCase.verifyEqual(actualB,expectedB);
        end

        function ordertest(testCase)
            V = vandermonde(2,2,[2 -1; 4 1],@grnlex);
            K = monomials(2,2,@grnlex);
            [actualA,actualB] = multmapnull(V,K,[2 1 0; 3 0 1],[1 2],...
                @grnlex);
            expectedA = [1 1; 2 4];
            expectedB = [1 11; 2 44];
            testCase.verifyEqual(actualA,expectedA);
            testCase.verifyEqual(actualB,expectedB);
        end

        function basistest(testCase)
            V = vandermonde(4,2,[2 -1]);
            K = monomials(4,2);
            [actualA,actualB] = multmapnull(V,K,[1 1 1],[5],...
                @chebyshev);
            expectedA = -2;
            expectedB = 2.5;
            testCase.verifyEqual(actualA,expectedA);
            testCase.verifyEqual(actualB,expectedB);
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
            V = vandermonde(5,3,[1 1 1; 2 2 2; 3 3 3]);
            [A1,B1] = multmapnull(V,K,ones(2,4),[1 2 4]);
            [A2,B2] = multmapnull(V,K,{ones(2,4)},[1 2 4]);
            [A3,B3] = multmapnull(V,K,{ones(2,4), ones(3,4)},[1 2 4]);
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
            Z = null(full(macaulay(system,3)));
            K = monomials(3,2);
            [A,B] = multmapnull(Z,K,[ones(2,1) eye(2)],[1 2 3 5]);
            actual = sort(eig(A,B));
            expected = [-1; -1; 1; 1];
            testCase.verifyLessThanOrEqual(actual-expected, ...
                testCase.tolerance);
        end

        function thesistest(testCase)
            p = [1 2 0; 1 1 1; -2 0 0];
            q = [1 0 2; 1 1 1; -2 0 0];
            system = systemstruct({p,q});
            Z = null(full(macaulay(system,4)));
            dgap = gap(Z,2);
            W11 = columncompr(Z,dgap,2);
            c = 1:3;
            [A,B] = multmapnull(W11,monomials(dgap,2),{[1 1 0],[1 0 1]},c);
            actual = sortrows([eig(pinv(A{1})*B{1}) eig(pinv(A{2})*B{2})]);
            expected = sortrows([1 1; -1 -1]);
            testCase.verifyLessThanOrEqual(actual-expected,testCase.tolerance);
        end
    end
end
