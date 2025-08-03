function [energy_total, unsafe, collision, profile, param] = v_sim(dt, beta, data, use_filter, vehicle_type, num_lead_vehicles, a_under)

    v1_ind = 6;
    a1_under = 8; % acceleration of the lead vehicle
    % a_under = 4; % acceleration of the ego vehicle
    % initiating: 
    % controller class: beta (vector, 1 by n, w_lead n by 1)
    % safety_filter class: a1_under, a_under (scalar)
    % car class: a_under (scalar)

    controller = v_controller(beta, vehicle_type); 
    safety_filter = v_safety_filter(a1_under,a_under); 
    car_module = v_dynamics(a_under, vehicle_type); 

    t_data = @(data,n) data.time{n};
    v_data = @(data,n) data.vel{n};
    h_data = @(data,n) data.hdwy{n};
    v_interp = @(t_data,v_data,time) interp1(t_data,v_data,time);
    
    % time for spectrum analysis
    time_data = t_data(data,v1_ind);
    time = time_data(1):dt:time_data(end);
    % V1 = v_interp(t_data(data,v1_ind), v_data(data,v1_ind), time);
    V_lead = zeros(num_lead_vehicles, length(time)); % [v1;v2;v3]
    lead_index = zeros(1,num_lead_vehicles);
    for i = v1_ind : -1: v1_ind - num_lead_vehicles + 1
        lead_index(v1_ind - i + 1) = i;
        V_lead(v1_ind - i + 1, :) = v_interp(t_data(data,i), v_data(data,i), time); % [v1;v2;v3]
    end
    
    % implemented sg filter to filter out the jittering in lead vehicle acceleration
    V1  = V_lead(1, :); % ego vehicle speed
    A1 = diff(V1)./dt;
    A1 = sgolayfilt(A1, 3, 21);

    % ego profile 
    h_ego = h_data(data, v1_ind + 1); 
    v_ego = v_data(data, v1_ind + 1);

    % IC
    H = time * 0;
    V = time * 0;
    Vh = time * 0;
    fuel_profile = time * 0;
    energy_profile = time * 0;

    ts = time(1:end-1); % for initiating the interval record
    a_ccc_profile = ts * 0; 
    a_cbf_profile = ts * 0;
    b_hat_profile = ts * 0;
    a_safe_profile = ts * 0;
    resist_profile = ts * 0;
    u_profile = ts * 0;
    u_ub_profile = ts * 0;
    u_lb_profile = ts * 0;
    dv_profile = ts * 0;

    % IC, using the initial starting point of the ego vehicle
    % v0 = V1(1);
    % h0 = controller.Vh_inverse(v0); % equilibrium headway
    H(1) = h_ego(1); % initial headway
    V(1) = v_ego(1); % initial speed
    fuel_total = 0;
    energy_total = 0;

    for i = 1:length(time)- 1
        v1 = V1(i); v = V(i); a1 = A1(i); h = H(i);
        v_lead = V_lead(:, i); % leading vehicle speed
        w_lead = min(v_lead, controller.v_max); 
        % controller -> safety filter -> powertrain constrain -> cost 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vh = controller.Vh(h);
        Vh(i) = vh;
        a_ccc = controller.ccc(vh, v, w_lead);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% safety filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [a_cbf, b_hat] = safety_filter.cbf(h, v, v1, a1);
        a_safe = min(a_ccc, a_cbf); % safety filter

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vehicle constrains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dh = v1 - v; 
        
        f_resist = car_module.resistance(v);
        [u_ub, u_lb] = car_module.bounds(v);
        
        if use_filter == 1
            u = a_safe + f_resist;
            u_sat = car_module.saturation(u, u_ub, u_lb);
            dv = -f_resist + u_sat;
        else
            u = a_ccc + f_resist;
            u_sat = car_module.saturation(u, u_ub, u_lb);
            dv = -f_resist + u_sat;
        end

        H(i+1) = h + dh * dt;
        V(i+1) = v + dv * dt;

        if V(i+1) < 0
            V(i+1) = 0; % avoid negative speed
        end
        if V(i+1) > controller.v_max
            V(i+1) = controller.v_max; % avoid over speed, range policy in controller
        end
        if i == length(time) - 1
            Vh(i+1) = controller.Vh(H(i+1));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cost %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [fuel, energy] = car_module.consumption(dv, u_sat, v, f_resist);
        fuel_total = fuel_total + fuel * dt;
        energy_total = energy_total + energy * dt;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% record %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fuel_profile(i+1) = fuel_total;
        energy_profile(i+1) = energy_total;
        a_ccc_profile(i) = a_ccc;
        a_cbf_profile(i) = a_cbf;
        b_hat_profile(i) = b_hat;
        a_safe_profile(i) = a_safe;
        resist_profile(i) = f_resist;
        u_profile(i) = u;
        u_ub_profile(i) = u_ub;
        u_lb_profile(i) = u_lb;
        dv_profile(i) = dv;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% safety and collision check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any((H - V .* safety_filter.tau) <0)
        unsafe = 1;
    else
        unsafe = 0;
    end

    if any(H <0)
        collision = 1;
    else
        collision = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % controller parameters
    param.beta = beta;
    param.alpha = controller.alpha;
    param.kappa = controller.kappa;
    param.h_st = controller.h_st;
    param.v_max = controller.v_max;
    % safety filter parameters
    param.tau = safety_filter.tau;
    param.gamma = safety_filter.gamma;  
    param.a1_under = a1_under;
    param.a_under = a_under;   
    % leading vehicle index
    param.lead_index = lead_index;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% profiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % states
    profile.H = H;
    profile.V = V;
    profile.V1 = V1;
    profile.Vh = Vh;
    profile.time = time;

    % feedback
    profile.A1 = A1;
    profile.dv = dv_profile;
    profile.u = u_profile;
    profile.u_ub = u_ub_profile;
    profile.u_lb = u_lb_profile;
    profile.a_ccc = a_ccc_profile;
    profile.a_cbf = a_cbf_profile;
    profile.b_hat = b_hat_profile;
    profile.a_safe = a_safe_profile;
    profile.resist = resist_profile;

    % cost
    profile.fuel = fuel_profile;
    profile.energy = energy_profile;

end
