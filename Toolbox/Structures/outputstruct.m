classdef outputstruct
    %OUTPUTSTRUCT   Class for gathering all the output details.
    %   obj = OUTPUTSTRUCT(time,nullity,shiftvalues) creates a structure 
    %   that contains the output information.
    %
    %   OUTPUTSTRUCT properties
    %       time         - necessary computation time.
    %       nullity      - nullity per degree.
    %       shiftvalues  - numerical values of the shift polynomial.
    %       multiplicity - flag to indicate multiplicities.
    %
    %   OUTPUTSTRUCT methods:
    %       disp - displays the output details.

    % Copyright (c) 2024 - Christof Vermeersch

    properties
        time 
        nullity
        shiftvalues
        multiplicity
    end

    methods
        function obj = outputstruct(time,nullity,shiftvalues)
            %OUTPUTSTRUCT   Constructor for the outputstruct class.
            %   obj = OUTPUTSTRUCT(time,nullity,shiftvalues) creates a 
            %   structure that contains the output information.
            
            % Store the different output properties:
            obj.time = time;
            obj.nullity = nullity;
            obj.shiftvalues = shiftvalues;
        end 

        function disp(~)
            %DISP   Information about the output.
            %   DISP(output) displays the output details.

            fprintf("     output structure \n")
        end
    end
end
