classdef testmacaulay < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function systemtest(testCase)
            sys = randsystem(2,2,2);
            M = full(macaulay(sys,2));
            N = [sys.coef{1}'; sys.coef{2}'];
            testCase.verifyLessThanOrEqual(M-N,testCase.tolerance);
        end

        function meptest(testCase)
            mep = randmep(2,2,3,2);
            M = full(macaulay(mep,2));
            N = [mep.mat{1} mep.mat{2} mep.mat{3} mep.mat{4} mep.mat{5} ...
                mep.mat{6}];
            testCase.verifyLessThanOrEqual(abs(M-N),testCase.tolerance);
        end

        function systemstructuretest(testCase)
            sys = systemstruct({[1 0 0; 2 0 1; 3 1 0; 4 2 0], ...
                [5 0 0; 6 0 1]});
            M = [1 2 3 0 0 4 0 0 0 0; 0 1 0 2 3 0 0 0 4 0; ...
                0 0 1 0 2 3 0 0 0 4; 5 6 0 0 0 0 0 0 0 0; ...
                0 5 0 6 0 0 0 0 0 0; 0 0 5 0 6 0 0 0 0 0; ...
                0 0 0 5 0 0 6 0 0 0; 0 0 0 0 5 0 0 6 0 0; ...
                0 0 0 0 0 5 0 0 6 0];
            N = macaulay(sys,3);
             testCase.verifyEqual(M,full(N));
        end

        function mepstructuretest(testCase)
            P = cell(6,1);
            for k = 1:6
                P{k} = randn(11,10);
            end
            Z = zeros(11,10);
            M = [P{1} P{2} P{3} P{4} P{5} P{6} Z Z Z Z; Z P{1} Z P{2} ...
                P{3} Z P{4} P{5} P{6} Z; Z Z P{1} Z P{2} P{3} Z P{4} ...
                P{5} P{6}];
            N = macaulay(mepstruct(P,[0 0; 0 1; 1 0; 0 2; 1 1; 2 0]),3);
            testCase.verifyLessThanOrEqual(abs(M-N),testCase.tolerance);
        end

        function complextest(testCase)
            sys = systemstruct({[1i 0 0; 2 0 1; 3-2i 1 0; 4 2 0], ...
                [5 0 0; 6 0 1]});
            M = [1i 2 3-2i 0 0 4 0 0 0 0; 0 1i 0 2 3-2i 0 0 0 4 0; ...
                0 0 1i 0 2 3-2i 0 0 0 4; 5 6 0 0 0 0 0 0 0 0; ...
                0 5 0 6 0 0 0 0 0 0; 0 0 5 0 6 0 0 0 0 0; ...
                0 0 0 5 0 0 6 0 0 0; 0 0 0 0 5 0 0 6 0 0; ...
                0 0 0 0 0 5 0 0 6 0];
            N = macaulay(sys,3);
             testCase.verifyEqual(M,full(N));
        end

        function ordertest(testCase)
            sys = systemstruct({[1 0 0; 2 1 0; 3 0 1; 4 0 2], ...
                [5 0 0; 6 1 0]});
            M = [1 2 3 0 0 4 0 0 0 0; 0 1 0 2 3 0 0 0 4 0; ...
                0 0 1 0 2 3 0 0 0 4; 5 6 0 0 0 0 0 0 0 0; ...
                0 5 0 6 0 0 0 0 0 0; 0 0 5 0 6 0 0 0 0 0; ...
                0 0 0 5 0 0 6 0 0 0; 0 0 0 0 5 0 0 6 0 0; ...
                0 0 0 0 0 5 0 0 6 0];
            N = macaulay(sys,3,@grnlex);
             testCase.verifyEqual(M,full(N));
        end

        function basistest(testCase)
            sys = systemstruct({[1 0 0; 2 0 1; 3 1 1; 4 2 0]});
            M = [1 2 0 0 3 4 0 0 0 0; 1 1 1.5 1 0 0 0 1.5 4 0; ...
                0 1.5 3 0 2 0 0 0 1.5 2];
            N = full(macaulay(sys,3,@chebyshev));
            testCase.verifyEqual(M,N);
        end

        function thesistest(testCase)
            p = [1 2 0; 1 0 2; -6 1 0; 7 0 0];
            q = [1 1 0; -1 0 1; -3 0 0];
            system = systemstruct({p,q});
            actual = macaulay(system,3,@grnlex);
            expected = [7 -6 0 1 0 1 0 0 0 0; 0 7 0 -6 0 0 1 0 1 0; ...
                0 0 7 0 -6 0 0 1 0 1; -3 1 -1 0 0 0 0 0 0 0; ...
                0 -3 0 1 -1 0 0 0 0 0; 0 0 -3 0 1 -1 0 0 0 0; ...
                0 0 0 -3 0 0 1 -1 0 0; 0 0 0 0 -3 0 0 1 -1 0; ...
                0 0 0 0 0 -3 0 0 1 -1];
            testCase.verifyEqual(full(actual),expected);
        end
    end
end
