%%% Paul Adkisson 
%%% 11.5.21
%%% Purpose Generate micro-stimulation current
function GenerateMicroStim(t, t_task, t_taskoff, stim_duration, stim_freq, ...
                          depol_block_thresh, depol_block_factor, pulse_amps, ...
                          dc_amps, N, num_group, brains, sim_path)
    I_ustim_base = zeros(length(t), N);
    dt = t(2) - t(1);
    stim_amps = [pulse_amps, dc_amps];
    lif_gl = 25*1e-9;
    lif_Cm = 0.5*1e-9;
    lif_tau = lif_Cm / lif_gl;
    C_300 = (1-exp(-stim_duration/lif_tau)); %time correction factor fro 300us/phase pulse
    EL = -70e-3; %mV
    Vs = -50e-3; %mV
    mir40 = mirror_est(40*1e-6);
    C_40 = (Vs - EL)/(C_300*mir40); %Threshold correction factor: 10uA --> AP at 40um
    
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        is_pulse = j <= length(pulse_amps);
        if is_pulse
            for i = 1:length(t) 
                if floor(1/(stim_freq*dt)) >= length(t)
                    mod_i = i;
                else
                    mod_i = mod(i, floor(1/(stim_freq*dt))) + 1;
                end
                if t(mod_i) <=  stim_duration
                    I_ustim_base(i, :) = stim_amp*C_40;
                elseif t(mod_i) > stim_duration && t(mod_i) <= 2*stim_duration
                    I_ustim_base(i, :) = -stim_amp*C_40;
                end
            end
        else
            I_ustim_base = ones(length(t), N)*stim_amp;
        end
        I_ustim_base(t<t_task|t>t_taskoff, :) = 0;
        for brain = brains
            brainpath = strcat(sim_path, sprintf("/brain%0.0f", brain));
            load(strcat(brainpath, "/r.mat"), "electric_r", "ball_r")
            I_ustim = [I_ustim_base(:, 1:num_group).*electric_r, zeros(length(t), N-num_group)];
            
            if is_pulse
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_pulse.mat", stim_amp*1e9)), "I_ustim")
            else
                task_time = t>=t_task & t<t_taskoff;
                super_thresh_neurons = I_ustim(t==t_task, :) > depol_block_thresh;
                super_thresh_amps = I_ustim(t==t_task, super_thresh_neurons);
                corrected_amps = super_thresh_amps - (super_thresh_amps - depol_block_thresh)*depol_block_factor;
                corrected_amps = repmat(corrected_amps, sum(task_time), 1);
                I_ustim(task_time, super_thresh_neurons) = corrected_amps; 
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_galvanic.mat", stim_amp*1e9)), "I_ustim")
            end
        end
        
        true_amps = I_ustim(t==t_task, 1:num_group);
        figure;
        hold on
        scatter(ball_r*1e6, true_amps*1e9)
        xlabel("Distance from Electrode (um)")
        ylabel("Stimulation Amplitude (nA)")
        if is_pulse
            title("Pulse")
        elseif stim_amp == 0
            title("Control")
        else
            title("Galvanic")
            yline(depol_block_thresh*1e9, 'r--')
            legend(["", "Depolarization Block Threshold"])
        end
        %}
    end
end

function Vm_ss = mirror_est(z)
    rho_e = 3; %extracellular resistivity (Ohm-m)
    N = 25; %Number of nodes
    n = (-(N-1)/2:(N-1)/2)'; %node numbers
    D = 10*10^(-6); %Fiber Diameter (m)
    delta_x = 100*D; %Internode Distance (m)
    x = delta_x*n; %node positions
    r = sqrt(x.^2 + z.^2); %distance from electrode to node (m)
    I_el = -10*1e-6;
    V_e = I_el*rho_e ./ (4*pi*z);
    V_es = I_el*rho_e ./ (4*pi*r);
    V_e_bar = mean(V_es);
    Vm_ss = (V_e_bar - V_e);
end