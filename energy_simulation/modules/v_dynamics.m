classdef v_dynamics
    properties
        % default value for heavy-duty truck
        m = 29484;
        k = 3.84;
        gamma = 0.006;
        pmax = 300.65e03;
        R = 0.504;
        I = 39.9;
        g = 9.81;
        a_max = 2;
        a  % lower bound acceleration
        % default hard coded values for fuel consumption model
        p2 = 1.8284;
        p1 = 0.0209;
        p0 = -0.1868;
        vehicle_type
    end
    methods
        function obj = v_dynamics(a, vehicle_type)
            obj.a = a;  
            obj.vehicle_type = vehicle_type;
        end
        %%%%%%%%%%%%%%%%%%%%%%% addubf more vehicle type for robustness %%%%%%%%%%%%%%%%%%%%%%%
        function f_resist = resistance(obj, v)
            if obj.vehicle_type == "car"
                f_resist = 0.0147 + 2.75e-04 * v^2;
            elseif obj.vehicle_type == "truck"
                meff = obj.m + obj.I / obj.R^2;
                b = obj.gamma * obj.m * obj.g / meff;
                k = obj.k / meff;
                f_resist = b + k * v^2;
            end
        end

        function [u_ub, u_lb] = bounds(obj,v)
            if obj.vehicle_type =="car"
                u_ub = min([obj.a_max, 0.285 * v + 2, -0.121 * v + 4.83]);
            elseif obj.vehicle_type == "truck"
                meff = obj.m + obj.I / obj.R^2;
                u_ub = min(obj.a_max, obj.pmax / (meff * v));
            end
            u_lb = -obj.a;
        end 

        %%%%%%%%%%%%%%%%%%%%%%%% for speed capping %%%%%%%%%%%%%%%%%%%%%%
        function u_sat = saturation(obj,u,u_ub,u_lb)
            u_sat = max(u_lb, min(u_ub, u));
        end

        %%%%%%%%%%%%%%%%%%%%%%%% for consumption %%%%%%%%%%%%%%%%%%%%%%
        function [fuel, energy] = consumption(obj, dv, u_sat, v, resist)
            % fuel part may not based on dv, confirm with the paper
            U1 = @(u, v) obj.p2 * u * v + obj.p1 * v + obj.p0;
            U2 = @(u, v) obj.p1 * v + obj.p0;
            if u_sat >= 0
                fuel = U1(u_sat,v);
            else
                fuel = U2(u_sat,v);
            end
            energy = v * max(0, dv + resist);
        end
    end
end
