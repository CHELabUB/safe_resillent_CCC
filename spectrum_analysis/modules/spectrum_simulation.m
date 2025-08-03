function [J, profile, param] = spectrum_simulation(data, dt, num_lead_vehicles, alpha, beta, kappa, sigma)
    % time for spectrum analysis
    v1_ind = 6; % index of the first vehicle in the data structure
    t_data = @(data,n) data.time{n};
    v_data = @(data,n) data.vel{n};
    v_interp = @(t_data,v_data,time) interp1(t_data,v_data,time);

    time_data = t_data(data,v1_ind);
    time = time_data(1):dt:time_data(end);
    V_lead = zeros(num_lead_vehicles, length(time)); % [v1;v2;v3]
    lead_index = zeros(1,num_lead_vehicles);
    for i = v1_ind : -1: v1_ind - num_lead_vehicles + 1
        lead_index(v1_ind - i + 1) = i;
        V_lead(v1_ind - i + 1, :) = v_interp(t_data(data,i), v_data(data,i), time); % [v1;v2;v3]
    end
    % pertubation
    V_avg = mean(mean(V_lead, 2));
    V_input = V_lead - V_avg;
    % from data
    magnitude_data = struct();
    phase_data = struct();
    for i  = 1:num_lead_vehicles
        dummy_name = sprintf('v%d', i);
        [magnitude_data.(dummy_name), phase_data.(dummy_name), omega] = fft_forward(dt, V_input(i, :)); 
    end

    % from Transfer function
    tf_magnitude_phase = struct();
    tf_spectrum = spectrum(alpha, beta, kappa, sigma); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spectrum reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % either sum^{1}{jj}( (sum^{1}{ii}((data_mag * TF_mag) * sin(data_phase + TF_phase)) )
    % or rewrite in paper format work (use a trigonometric identity)
    V_recon = zeros(size(time));
    magnitude_tf = zeros(num_lead_vehicles, length(omega));
    phase_tf = zeros(num_lead_vehicles, length(omega));
    magnitude_recon = zeros(1, length(omega));
    phase_recon = zeros(1, length(omega));
    J = 0;
    for ii = 1:length(omega)
        w = 2 * pi * omega(ii);
        magnitude_hat  = zeros(1, num_lead_vehicles); 
        phase_hat = zeros(1, num_lead_vehicles); 
        for jj = 1:num_lead_vehicles
            dummy_name = sprintf('TF%d', jj);
            % Calculate the transfer function magnitude and phase for each vehicle
            [tf_magnitude_phase.(dummy_name).magnitude, tf_magnitude_phase.(dummy_name).phase] = ...
                tf_spectrum.TF_magnitude_phase(w * 1i, jj);
            magnitude_hat(jj) = magnitude_data.(sprintf('v%d', jj))(ii) * tf_magnitude_phase.(dummy_name).magnitude;
            phase_hat(jj) = phase_data.(sprintf('v%d', jj))(ii) + tf_magnitude_phase.(dummy_name).phase;
            magnitude_tf(jj, ii) = tf_magnitude_phase.(dummy_name).magnitude;
            phase_tf(jj, ii) = tf_magnitude_phase.(dummy_name).phase;
        end

        % reconstruction
        magnitude_recon(ii) = sqrt((magnitude_hat * cos(phase_hat'))^2 + ...
                            (magnitude_hat * sin(phase_hat'))^2); % Magnitude of the reconstructed signal
        phase_recon(ii) = atan2(magnitude_hat * sin(phase_hat'), magnitude_hat * cos(phase_hat')); 

        % signal reconstruct
        V_recon = V_recon + magnitude_recon(ii) * sin(w * (time - time(1)) + phase_recon(ii));
        % Energy computing
        J = J + w^2 * magnitude_recon(ii)^2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% record %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % params
    param.alpha = alpha;
    param.beta = beta;
    param.kappa = kappa;
    param.sigma = sigma;
    param.num_lead_vehicles = num_lead_vehicles;
    param.lead_index = lead_index;

    % profiles
    % data (magnitude, phase)
    profile.magnitude_data = magnitude_data;
    profile.phase_data = phase_data;
    % transfer function (magnitude, phase)
    profile.magnitude_tf = magnitude_tf;
    profile.phase_tf = phase_tf;
    % signal hat (magnitude, phase)
    profile.magnitude_hat = magnitude_hat;
    profile.phase_hat = phase_hat;
    % signal for reconstruction (magnitude, phase)
    profile.magnitude_recon = magnitude_recon;
    profile.phase_recon = phase_recon;
    % reconstructed signal 
    profile.V_lead = V_lead;
    profile.V_recon = V_recon + V_avg;
    profile.time = time;
    
end
