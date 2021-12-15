clear;
sim_name = "Test";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
default_colors = get(gca, "colororder");

%{
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
%}

% Spatial Electrode Effects
min_x = sqrt(50*1e-12); %minimum distance of 10um
max_x = 2*1e-3; %From Levitt et al.
min_y = sqrt(50*1e-12);
max_y = 200*1e-6; %From Levitt et al.
regular_x = (max_x - min_x) / num_group;
regular_y = (max_y - min_y) / num_group;
group_idx = 0:num_group-1;
for brain = brains
    rpath = strcat(sim_path, sprintf("/brain%0.0f/r.mat", brain));
    load(rpath, "ball_r")
    rng(brain);
    ball_y = min_y + (min(ball_r, max_y) - min_y).*rand(1, num_group);
    ball_x = sqrt(ball_r.^2 - ball_y.^2);
    neg_y = rand(1, num_group) <= 0.5;
    neg_x = rand(1, num_group) <= 0.5;
    ball_y(neg_y) = -1*ball_y(neg_y);
    ball_x(neg_x) = -1*ball_x(neg_x);
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j<=length(pulse_amps);
        if pulse
            decpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse/decisions.mat", ...
                [sim_name, brain, stim_amp*1e9]);
        else
            decpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic/decisions.mat", ...
                [sim_name, brain, stim_amp*1e9]);
        end
        load(decpath, 'totspikes_g1')
        figure;
        hold on
        scatter(ball_x, ball_y, [], totspikes_g1)
        scatter(0, 0, "k", "filled")
        hold off
        xlim([-max_x, max_x])
        ylim([-max_x, max_x])
        colormap jet
        colorbar
        if pulse
            title(sprintf("%0.1fnA Pulse", ...
                [stim_amp*1e9]))
        else
            title(sprintf("%0.1fnA GVS, Coherence: %0.1f%%", ...
                [stim_amp*1e9]))
        end
    end
end
