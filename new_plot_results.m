clear;
sim_name = "Test";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
default_colors = get(gca, "colororder");

%Example FRs
ex_c = -100/100;
ex_trial = 1;
ex_brain = 1;
figure;
axs = zeros(length(stim_amps), 1);
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    if pulse
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
            [sim_name, ex_brain, stim_amp*1e9]);
    else
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
            [sim_name, ex_brain, stim_amp*1e9]);
    end
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), "pop_frs")
    axs(j) = subplot(length(stim_amps), 1, j);
    hold on
    for i = 1:p+2
        plot(t, pop_frs(:, i))
    end
    hold off
    if pulse
        title(sprintf("%0.1fnA Pulse, Coherence: %0.1f%%", ...
            [stim_amp*1e9, ex_c*100]))
    else
        title(sprintf("%0.1fnA GVS, Coherence: %0.1f%%", ...
            [stim_amp*1e9, ex_c*100]))
    end
    if j == length(stim_amps)
        xlabel("Time (s)")
    elseif j == ceil(length(stim_amps)/2)
        ylabel("Population Firing Rate (Hz)")
    elseif j == 1
        legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory"], "Location", "northwest")
    end
end
linkaxes(axs)