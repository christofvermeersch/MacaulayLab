classdef solutionstruct
    %SOLUTIONSTRUCT   Class for gathering all the solution data.
    %   obj = SOLUTIONSTRUCT(values) creates a structure that contains
    %   all the information about the solution(s).
    %
    %   SOLUTIONSTRUCT properties:
    %       values - numerical values.
    %
    %   SOLUTIONSTRUCT methods:
    %       num  - return the numerical values of the solution(s).
    %       disp - displays the solution(s) the problem.
    %       info - gives a more elaborate description of the solution(s).

    % Copyright (c) 2024 - Christof Vermeersch

    properties
        values
        ma
        mb
        res
        maxres
    end

    methods
        function obj = solutionstruct(values)
            %SOLUTIONSTRUCT   Constructor for the solutionstruct class.
            %   obj = SOLUTIONSTRUCT(values) creates a structure that 
            %   contains all the information about the solutions.

            obj.values = values;
        end 

        function values = num(obj)
            %NUM   Numerical values of the solution(s).
            %   values = NUM(solutions) returns the numerical values of the
            %   solution(s).

            values = obj.values;
        end

        
        function disp(obj)
            %DISP   Solution(s) of a problem.
            %   DISP(solutions) displays the solution(s) the problem.

            disp(num(obj))
        end

        function output = info(obj)
            %INFO   Information about the solution(s).
            %   INFO(solutions) gives a more elaborate description of the 
            %   solution(s).
            
            text = sprintf("     total number of solutions: %d \n", ...
                obj.mb);
            text = text + ...
                sprintf("     number of affine solutions: %d \n", obj.ma);

            if nargout > 0
                output = text;
            else
                fprintf(text);
            end
        end
    end
end