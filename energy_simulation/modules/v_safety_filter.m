classdef v_safety_filter
    properties
        a1 
        a 
        tau = 1;
        gamma = 1.8;
    end
    methods

        function obj = v_safety_filter(a1, a)
            obj.a1 = a1; 
            obj.a = a;   
        end

        function [a_cbf,barrier] = cbf(obj, h, v, v1, a1)
            B1 = @(v1,v) v * obj.tau;
            B2 = @(v1,v) v * obj.tau + (v - obj.a * obj.tau)^2 / (2 * obj.a)...
                - v1^2 / (2 * obj.a1);
            f1 = @(v) sqrt(obj.a1 / obj.a) * (v - obj.a * obj.tau);
            dB1dv = @(v1,v) obj.tau;
            dB1dv1 = @(v1,v) 0;
            dB2dv = @(v1,v) v / obj.a;
            dB2dv1 = @(v1,v) -v1 / obj.a1;

            % safety filter 
            a_barrier1 = @(v1,v,h,a1) 1 / dB1dv(v1,v) * (v1 - v - dB1dv1(v1,v) * a1 + ...
                obj.gamma * (h - B1(v1,v)));
            a_barrier2 = @(v1,v,h,a1) 1 / dB2dv(v1,v) * (v1 - v - dB2dv1(v1,v) * a1 + ...
                obj.gamma * (h - B2(v1,v)));

            if f1(v) > v1
                a_cbf = a_barrier2(v1,v,h,a1);
                barrier = B2(v1,v);
            else
                a_cbf = a_barrier1(v1,v,h,a1);
                barrier = B1(v1,v);
            end
        end
    end
end
