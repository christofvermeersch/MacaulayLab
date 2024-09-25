classdef testmultiplicity < matlab.unittest.TestCase
    methods (Test)
        % Test the method:
        function multiplicity1test(testCase)
            X = [0.9999; 1.0001; 2];
            C = multiplicity(X);
            testCase.verifyEqual(C,[1; 1; 2]);
        end

        function multiplicity2test(testCase)
            X = [0.9999 1.0001; 1.0001 0.9999; 2 3];
            C = multiplicity(X);
            testCase.verifyEqual(C,[1 1; 1 1; 2 3]);
        end

        function complexitytest(testCase)
            X = [0.9999+1i 1.0001-1i; 1.0001+1i 0.9999-1i; 2 3];
            C = multiplicity(X);
            testCase.verifyEqual(C,[1+1i 1-1i; 1+1i 1-1i; 2 3]);
        end

        function clustertest(testCase)
            X = randn(5,3);
            X(2,:) = X(1,:);
            X(5,:) = X(3,:)*1.0001;
            [~,E] = multiplicity(X);
            testCase.verifyEqual(E(1),E(2));
            testCase.verifyEqual(E(3),E(5));
        end
    end
end
