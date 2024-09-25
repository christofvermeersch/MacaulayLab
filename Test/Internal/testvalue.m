classdef testvalue < matlab.unittest.TestCase
    properties (Constant)
        tolerance = 1e-10;
    end

    methods (Test)
        % Test the method:
        function univariatetest(testCase)
            p = [1 3; 2 2; 3 1; 4 0];
            system = systemstruct({p});
            x = -2;
            expected = -2;
            actual = value(x,system);
            testCase.verifyEqual(actual,expected);
        end

        function complextest(testCase)
            p = [1 3; 2i 2; 3-2i 1; 4 0];
            system = systemstruct({p});
            x = -1-2i;
            expected = -8i;
            actual = value(x,system);
            testCase.verifyEqual(actual,expected);
        end

        function bivariatetest(testCase)
            p = [1 1 1; 2 2 0; 3 2 2];
            system = systemstruct({p});
            x = [1 2];
            expected = 16;
            actual = value(x,system);
            testCase.verifyEqual(actual,expected);
        end


        function polynomialtest(testCase)
            p = {[1 0 0 0; 2 1 0 0; 3 0 1 1; 4 0 3 0; 5 2 2 2], ...
                [1 0 0; 1 1 0; 1 0 1; 1 2 0; 1 1 1; 1 0 2] ...
                [1 0; 2 1; 1 2]};
            system = systemstruct(p);
            x = {[1 -2 3], [0.5 2], 2};

            expected1 = {133, 8.75, 9};
            expected2 = {476, 11, 12};            

            for k = 1:3
                val = value(x{k},system.coef{k},system.supp{k});
                testCase.verifyLessThanOrEqual(val-expected1{k}, ...
                    testCase.tolerance);
            end

            for k = 1:3
                val = value(x{k},system.coef{k},system.supp{k},@chebyshev);
                testCase.verifyLessThanOrEqual(val-expected2{k}, ...
                    testCase.tolerance);
            end
        end 

        function meptest(testCase)
            P = {{[1 2; 3 4; 3 4]; [2 1; 0 1; 1 3]; zeros(3,2); ...
                zeros(3,2); [3 4; 2 1; 0 1]; [1 2; 4 2; 2 1];}, ...
                {ones(5,5); ones(5,5); ones(5,5); ones(5,5)}};
            n = {2, 3};
            dmax = {2, 1};
            x = {[2 2], [1 -2 3]};

            expected1 = {[21 28; 27 18; 13 18], 3*ones(5,5)};
            expected2 = {[24 34; 39 24; 19 21], 3*ones(5,5)};            

            for k = 1:2
                mep = mepstruct(P{k},dmax{k},n{k});
                val = value(x{k},mep);
                testCase.verifyLessThanOrEqual(val-expected1{k}, ...
                    testCase.tolerance);
            end

            for k = 1:2
                mep = mepstruct(P{k},dmax{k},n{k});
                val = value(x{k},mep,@chebyshev);
                testCase.verifyLessThanOrEqual(val-expected2{k}, ...
                    testCase.tolerance);
            end
        end
    end
end
