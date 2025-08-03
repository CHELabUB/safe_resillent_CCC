classdef spectrum
    properties
        alpha
        beta
        kappa 
        delay
    end
    
    methods
        function obj = spectrum(alpha, beta, kappa, delay)
            obj.alpha = alpha; 
            obj.beta = beta; 
            obj.kappa = kappa; 
            obj.delay = delay; 
        end
        
        function G = TF(obj, w, num_vehicle)
            % Transfer function G(jw)
            if num_vehicle == 1
                G = (obj.alpha * obj.kappa + w * obj.beta(num_vehicle)) / ...
                (exp(obj.delay * w) * w^2 + (obj.alpha + sum(obj.beta)) * w + obj.alpha * obj.kappa);
            else 
                G  = obj.beta(num_vehicle) * w / ...
                (exp(obj.delay * w) * w^2 + (obj.alpha + sum(obj.beta)) * w + obj.alpha * obj.kappa);
            end
        end

        function [G_magnitude, G_phase] = TF_magnitude_phase(obj, w, num_vehicle)
            % Calculate magnitude and phase of the transfer function
            G = obj.TF(w, num_vehicle);
            G_magnitude = abs(G);
            G_phase = angle(G); % Phase in radians
        end
    end
end
