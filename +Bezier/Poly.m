classdef Poly
    methods (Static)
        function V = hyp2vert(A,b)
            Vert = cddmex('extreme',struct('A',A,'B',b));
            V = Vert.V;
        end

        function [A,b] = vert2hyp(V)
            H = cddmex('hull',struct('V',V));
            A = H.A;
            b = H.B;
        end

        function V = conv(V)
            if ~isempty(V)
                ind = convhull(V);
                V = V(ind,:);
            end
        end
    end
end
