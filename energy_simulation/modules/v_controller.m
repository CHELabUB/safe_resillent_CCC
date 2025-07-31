classdef v_controller
    properties
        % for CCC controller
        alpha = 0.4;
        beta 
        % for range policy, Vh result also for CCC controller
        kappa = 0.6;
        h_st = 5;
        v_max = 30;
    end

    methods

        function obj = v_controller(beta, vehicle_type)
            obj.beta = beta; 
            if vehicle_type == "car"
                obj.v_max = 35;
            end
        end

        function vh = Vh(obj,h)
            h_go = 1 / obj.kappa * obj.v_max + obj.h_st;
            if h <=  obj.h_st
                vh = 0;
            elseif h > obj.h_st && h < h_go
                vh = obj.kappa * (h - obj.h_st);
            else
                vh = obj.v_max;
            end
        end

        function h_inverse = Vh_inverse(obj,v)
            % possible use: when v is given, back calculate h
            % h_inverse = 1 / obj.kappa * obj.v_max + obj.h_st;
            h_inverse = 1 / obj.kappa * v + obj.h_st;
        end

        function a_ccc = ccc(obj, vh, v, v_lead)
            % Note: 
            % vh is scaler, v is scaler, a_ccc is scaler
            % size(beta) = (1,3), while size (v_lead - v) = (3,1)
            a_ccc = obj.alpha * (vh - v) + obj.beta * (v_lead - v);
        end
    end
end
