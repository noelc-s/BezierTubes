classdef LinearSystem

    properties
        A;
        B;
    end

    methods
        function obj = LinearSystem(A,B)
            obj.A = A;
            obj.B = B;
        end
    end

end