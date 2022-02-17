%%% Paul Adkisson
%%% 2/14/2022
%%% Plot Synchrony
function plot_sync(pulse_amps, stim_amps, t, t_task, t_taskoff, stim_freq, num_group, ...
                        default_colors, brains, num_brains, N_start, N_end, ...
                        win_start, win_stop, c_win, c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials)
    for sim_name = ["EMBC Disconnected", "EMBC I_b100"]
        disp(sim_name)

        stim_sync = zeros(length(stim_amps), num_brains, num_group, num_group);
        for brain = brains
            fprintf("brain %0.0f \n", brain)
            for j = 1:length(stim_amps)
                stim_amp = stim_amps(j);
                pulse = j<=length(pulse_amps);
                if pulse
                    output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                        [sim_name, brain, stim_amp*1e9]);
                    stim_coherences = pulse_coherences;
                    disp("Pulsatile")
                else
                    output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                        [sim_name, brain, stim_amp*1e9]);
                    if stim_amp == 0
                        stim_coherences = control_coherences;
                        if brain ~= 1
                            continue
                        end
                        disp("Control")
                    else
                        stim_coherences = galvanic_coherences;
                        disp("Galvanic")
                    end
                end
                %load(strcat(output_stimpath, "/decisions.mat"), "final_decisions")
                load(sprintf("Simulation %s/brain%0.0f/data/0.0nA_galvanic/decisions.mat", ...
                        [sim_name, brain]), 'final_decisions'); %using only control final decs
                if sim_name == "EMBC I_b100"
                    num_wins = sum(final_decisions(start_trial:end_trial, stim_coherences==0)==1, 'all') * ones(num_group, num_group);
                else
                    num_wins = num_trials * ones(num_group, num_group);
                end
                nan_neurons = zeros(num_trials, num_group);
                for trial = start_trial:end_trial
                    fprintf("Trial %0.0f \n", trial)
                    if sim_name=="EMBC I_b100" && final_decisions(trial, stim_coherences==c) ~= 1
                        disp("P1 lost")
                        continue %skip trials where P1 doesn't win
                    end
                    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [c, trial])), ...
                        "recspikes")
                    %pairwise_sync = get_pairwise_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
                    pairwise_sync = get_sym_sync(recspikes, N_start, N_end, t, c_win, win_start, win_stop);
                    no_spike_neurons = all(isnan(pairwise_sync), 1);
                    num_wins(no_spike_neurons, :) = num_wins(no_spike_neurons, :) - 1;
                    num_wins(:, no_spike_neurons) = num_wins(:, no_spike_neurons) - 1;
                    pairwise_sync = pairwise_sync(~no_spike_neurons, ~no_spike_neurons);
                    stim_sync(j, brain, ~no_spike_neurons, ~no_spike_neurons) = ...
                        reshape(stim_sync(j, brain, ~no_spike_neurons, ~no_spike_neurons), size(pairwise_sync)) + pairwise_sync;
                    nan_neurons(trial, :) = no_spike_neurons;
                end
                if any(all(nan_neurons==1, 1)) %no spikes for all 36 trials
                    stim_sync(j, brain, all(nan_neurons==1, 1), :) = NaN;
                    stim_sync(j, brain, :, all(nan_neurons==1, 1)) = NaN;
                end
                stim_sync(j, brain, :, :) = reshape(stim_sync(j, brain, :, :), size(num_wins)) ./ num_wins;
            end
        end

        pulse_sync = reshape(mean(stim_sync(1, :, :, :), 2, 'omitnan'), [num_group, num_group]);
        galvanic_sync = reshape(mean(stim_sync(2, :, :, :), 2, 'omitnan'), [num_group, num_group]);
        control_sync = reshape(stim_sync(3, 1, :, :), [num_group, num_group]);

        if sim_name == "EMBC I_b100"
            save("matdata/connected.mat", "pulse_sync", "galvanic_sync", "control_sync")
        else
            save("matdata/disconnected.mat", "pulse_sync", "galvanic_sync", "control_sync")
        end

        nan_color = uint8([0, 0, 128]);
        ticks = log10([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]);
        tick_labels = ["1", "10", "", "", "", "", "", "", "", "", "100"];
        figure;
        h = heatmap(log10(pulse_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Pulsatile")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        figure;
        h = heatmap(log10(galvanic_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Galvanic")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        figure;
        h = heatmap(log10(control_sync*100), 'ColorLimits', [0, 2], 'MissingDataColor', nan_color);
        colormap(hot)
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Control")
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = ticks;
        cb.TickLabels = tick_labels;

        %Calculate fraction of neuron pairs that are "significantly
        %synchronized" (>3SD from control mean)
        ctrl_mean = mean(control_sync, 'all', 'omitnan')
        ctrl_std = std(control_sync, [], 'all', 'omitnan')
        ctrl_synced = double(control_sync > ctrl_mean + 2.5*ctrl_std);
        ctrl_perc_synced = sum(ctrl_synced, 'all') / (length(ctrl_synced)^2 - num_group)
        ps_perc_synced = zeros(num_brains, 1); 
        gs_brainsync = zeros(num_brains, 1);
        for brain = brains
            ps_synced = stim_sync(1, brain, :, :) > ctrl_mean + 2.5*ctrl_std;
            ps_perc_synced(brain) = sum(ps_synced, 'all') / (length(ps_synced)^2 - num_group);
            gs_synced = stim_sync(2, brain, :, :) > ctrl_mean + 2.5*ctrl_std;
            gs_perc_synced(brain) = sum(gs_synced, 'all') / (length(gs_synced)^2 - num_group);
        end
        norm_ps_synced = ps_perc_synced - ctrl_perc_synced;
        norm_gs_synced = gs_perc_synced - ctrl_perc_synced;
        [~, p_val] = ttest2(norm_ps_synced, norm_gs_synced)
        mean_norm_ps_sync = mean(norm_ps_synced)
        sem_norm_ps_sync = std(norm_ps_synced) / sqrt(num_brains)
        mean_norm_gs_sync = mean(norm_gs_synced)
        sem_norm_gs_sync = std(norm_gs_synced) / sqrt(num_brains)

        %average for plotting
        ps_synced = double(pulse_sync > ctrl_mean + 2.5*ctrl_std);
        gs_synced = double(galvanic_sync > ctrl_mean + 2.5*ctrl_std);
        figure;
        h = heatmap(ps_synced);
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        colormap(hot)
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Pulse")

        figure;
        h = heatmap(gs_synced);
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        colormap(hot)
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Galvanic")

        figure;
        h = heatmap(ctrl_synced);
        h.XDisplayLabels = nan(size(h.XDisplayData));
        h.YDisplayLabels = nan(size(h.YDisplayData));
        colormap(hot)
        xlabel("neuron 2")
        ylabel("neuron 1")
        title("Control")
    end
end