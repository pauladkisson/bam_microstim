%%% Paul Adkisson 
%%% 11.5.21
%%% Purpose Generate micro-stimulation current
function GenerateMicroStim(t, t_task, t_taskoff, stim_duration, stim_freq, ...
                           pulse_amps, dc_amps, N, num_group, num_brains, ...
                           sim_path)
    I_ustim_base = zeros(length(t), N);
    dt = t(2) - t(1);
    stim_amps = [pulse_amps, dc_amps];
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
                    I_ustim_base(i, :) = stim_amp;
                elseif t(mod_i) > stim_duration && t(mod_i) <= 2*stim_duration
                    I_ustim_base(i, :) = -stim_amp;
                end
            end
        else
            I_ustim_base = ones(length(t), N)*stim_amp;
        end
        I_ustim_base(t<t_task|t>t_taskoff, :) = 0;
        for brain = 1:num_brains
            brainpath = strcat(sim_path, sprintf("/brain%0.0f", brain));
            load(strcat(brainpath, "/r.mat"), "electric_r")
            I_ustim = [I_ustim_base(:, 1:num_group).*electric_r, zeros(length(t), N-num_group)];
            
            if is_pulse
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_pulse.mat", stim_amp*1e9)), "I_ustim")
            else
                basepath = strcat(brainpath, "/ustim");
                mkdir(basepath)
                save(strcat(basepath, sprintf("/%0.1fnA_galvanic.mat", stim_amp*1e9)), "I_ustim")
            end
        end

        %{
        ustim_desired = I_ustim(t==t_task, 1:num_group)*1e9;
        mean_desired = mean(ustim_desired);
        one_pA_percent_desired = sum(ustim_desired>=1) / num_group;
        fprintf("Mean Stimulation: %0.1f nA \n", mean_desired)
        fprintf("%% of Neurons > 1pA: %0.1f \n", one_pA_percent_desired*100)
        
        ustim_undesired = I_ustim(t==t_task, num_group+1:end)*1e9;
        mean_undesired = mean(ustim_undesired);
        one_pA_percent_undesired = sum(ustim_undesired>=1) / (N-num_group);
        fprintf("Mean Stimulation: %0.1f nA \n", mean_undesired)
        fprintf("%% of Neurons > 1pA: %0.1f \n", one_pA_percent_undesired*100)
        
        figure;
        plot(t, I_ustim(:, 1)*1e12)
        title(sprintf("Stim Amp: %0.0fnA", stim_amp*1e9))
        
        figure;
        hold on
        scatter(electric_r(1:num_group), ustim_desired)
        scatter(electric_r(num_group+1:end), ustim_undesired)
        hold off
        legend(["Desired", "Undesired"]) 
        %}
    end
end