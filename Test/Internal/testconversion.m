classdef testconversion < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function conversiontest(testCase)
            syms a g1 g2

            p(1) = 4*a - 2*a*g1 - a^2*g2 - a^3*g1 + 3;
            p(2) = 3*a - g1 - 2*a*g2 - 2*a^2*g1 - a^3*g2 + 2;
            p(3) = 2*a - g2 - a*g1 - 2*a^2*g2 + 1;

            system = systemstruct(p);
            mep = conversion(system,[2,3]);
            
            M = cell(4,1);
            M{1} = [3 0 0; 2 -1 0; 1 0 -1];
            M{2} = [4 -2 0; 3 0 -2; 2 -1 0];
            M{3} = [0 0 -1; 0 -2 0; 0 0 -2];
            M{4} = [0 -1 0; 0 0 -1; 0 0 0];
            m = false(4,1);

            N = cell(4,1);
            N{1} = [3 0 0; 2 0 -1; 1 -1 0];
            N{2} = [4 0 -2; 3 -2 0; 2 0 -1];
            N{3} = [0 -1 0; 0 0 -2; 0 -2 0];
            N{4} = [0 0 -1; 0 -1 0; 0 0 0];
            n = false(4,1);

            for i = 1:4
                A = mep.mat{i};
                B = M{i};
                C = N{i};
                m(i) = all(A(:) == B(:));
                n(i) = all(A(:) == C(:));
            end

            testCase.verifyTrue(all(or(m,n)))
        end
    end
end
