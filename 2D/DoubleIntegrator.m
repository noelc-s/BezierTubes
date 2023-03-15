classdef DoubleIntegrator


    properties
        A;
        B;
    end

methods
    function obj = DoubleIntegrator()
        A = [0 1; 0 0];
        B = [0; 1];
        obj.A = A;
        obj.B = B;
    end
end

end