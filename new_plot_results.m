clear;
sim_name = "EMBC";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
default_colors = get(gca, "colororder");
num_trials = end_trial - start_trial + 1;
brains = 1:2;

%{
%Example FRs
ex_c = -3.2/100;
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

%{
%Example Raster: P1 vs. Distance & all groups vs. total activity (ranked) 
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = -0.512;
ex_brain = 1;
load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, ex_brain]), "ball_r")
p1_only = figure;
all_groups = figure;
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
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
        "recspikes")
    spikes = zeros(length(t), N);
    neuron_num = 1:N;
    for nn = 1:N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    
    g1_spikes = spikes(:, 1:num_group);
    [g1_time_idx, g1_neuron_idx, g1_idx] = get_spike_idx(g1_spikes);
    g2_spikes = spikes(:, num_group+1:2*num_group);
    [g2_time_idx, g2_neuron_idx, g2_idx] = get_spike_idx(g2_spikes);
    ns_spikes = spikes(:, 2*num_group+1:N_E);
    [ns_time_idx, ns_neuron_idx, ns_idx] = get_spike_idx(ns_spikes);
    int_spikes = spikes(:, N_E+1:end);
    [int_time_idx, int_neuron_idx, int_idx] = get_spike_idx(int_spikes);
    
    figure(p1_only);
    subplot(length(stim_amps), 1, j)
    scatter(t(g1_time_idx), ball_r(g1_idx(g1_neuron_idx)), "Marker", "|", ...
        "MarkerFaceColor", default_colors(1, :), "MarkerEdgeColor", default_colors(1, :))
    xlabel("Time (s)")
    ylabel("Distance from Electrode (um)")
    if pulse
        title(sprintf("Pulse, c=%0.1f", ex_c*100))
    else
        if abs(stim_amp) > 0
            title(sprintf("Galvanic, c=%0.1f", ex_c*100))
        else
            title(sprintf("Control, c=%0.1f", ex_c*100))
        end
    end
    
    figure(all_groups)
    subplot(length(stim_amps), 4, (j-1)*4 + 1)
    scatter(t(g1_time_idx), g1_neuron_idx, "Marker", "|", ...
        "MarkerFaceColor", default_colors(1, :), "MarkerEdgeColor", default_colors(1, :))
    if j == 1
        ylabel([""; "Pulse"])
        xlabel("")
        title([""; "Population 1"])
    elseif j == 2
        ylabel(["Neuron Number (ranked by total activity)"; "Galvanic"])
    else
        ylabel([""; "Control"])
    end
    subplot(length(stim_amps), 4, (j-1)*4 + 2)
    scatter(t(g2_time_idx), (g2_neuron_idx), "Marker", "|", ...
        "MarkerFaceColor", default_colors(2, :), "MarkerEdgeColor", default_colors(2, :))
    if j == 1
        title([sprintf("Coherence %0.1f%%", ex_c*100); "Population 2"])
    end
    subplot(length(stim_amps), 4, (j-1)*4 + 3)
    scatter(t(ns_time_idx), (ns_neuron_idx), "Marker", "|", ...
        "MarkerFaceColor", default_colors(3, :), "MarkerEdgeColor", default_colors(3, :))
    if j == length(stim_amps)
        xlabel("Time (s)")
    elseif j == 1
        title([""; "Non-Selective"])
    end
    subplot(length(stim_amps), 4, (j-1)*4 + 4)
    scatter(t(int_time_idx), (int_neuron_idx), "Marker", "|", ...
        "MarkerFaceColor", default_colors(4, :), "MarkerEdgeColor", default_colors(4, :))
    if j == 1
        title([""; "Interneurons"])
    end
end
%}

%{
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
%}


%Accuracies
c = 0:0.01:1;
figure;
hAx = axes;
%hAx.XScale = 'log';
hold on
idx = 1;
fig_leg = [];
for brain = brains
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
            stim_coherences = pulse_coherences;
            if brain == 1
                fig_leg = [fig_leg, sprintf("%0.1fnA Pulse", stim_amp*1e9)];
            end
        else
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
            if stim_amp == 0
                stim_coherences = control_coherences;
            else
                stim_coherences = galvanic_coherences;
            end
            if brain == 1
                fig_leg = [fig_leg, sprintf("%0.1fnA Galvanic", stim_amp*1e9)];
            end
        end
        load(strcat(datapath, "/decisions.mat"), "avg_acc", "decisions");
        %coeffs
        %scatter(coherences, avg_acc, dc_colors(idx))
        %plot(c, weibull(coeffs, c), dc_colors(idx))
        plot(coherences, avg_acc, 'o-', 'Color', default_colors(j, :))
    end
end
hold off
xlabel("Coherence")
ylabel("Accuracy")
legend(fig_leg)
%f = flipud(get(gca, 'Children'));
%legend([f(2), f(4), f(6)], "I-dc=-4pA", "I-dc=0pA", "I-dc=4pA")
%legend(compose("I-dc=%0.0fpA", I_dcs(1, :)*1e12))
%}

%Final Accuracies
c = 0:0.01:1;
figure;
hAx = axes;
%hAx.XScale = 'log';
hold on
idx = 1;
fig_leg = [];
for brain = brains
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
            if brain == 1
                fig_leg = [fig_leg, sprintf("%0.1fnA Pulse", stim_amp*1e9)];
            end
        else
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
            if brain == 1
                fig_leg = [fig_leg, sprintf("%0.1fnA Galvanic", stim_amp*1e9)];
            end
        end
        load(strcat(datapath, "/decisions.mat"), "avg_final_acc");
        %coeffs
        %scatter(coherences, avg_acc, dc_colors(idx))
        %plot(c, weibull(coeffs, c), dc_colors(idx))
        plot(coherences, avg_final_acc, 'o-', 'Color', default_colors(j, :))
    end
end
hold off
xlabel("Coherence")
ylabel("Final Accuracy")
legend(fig_leg)
%}

%{
%Decision Times
figure;
hAx = axes;
%hAx.XScale = 'log';
hold on
idx = 1;
fig_leg = [];
for brain = brains
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
            if brain == 1
                fig_leg = [fig_leg, sprintf("%0.1fnA Pulse", stim_amp*1e9)];
            end
        else
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
            if brain == 1
                fig_leg = [fig_leg, sprintf("%0.1fnA Galvanic", stim_amp*1e9)];
            end
        end
        load(strcat(datapath, "/decisions.mat"), "avg_dts", "std_dts");
        %coeffs
        %scatter(coherences, avg_acc, dc_colors(idx))
        %plot(c, weibull(coeffs, c), dc_colors(idx))
        errorbar(coherences, avg_dts, std_dts./sqrt(num_trials), 'o-', ...
            'Color', default_colors(j, :))
    end
end
hold off
xlabel("Coherence")
ylabel("Decision Time (s)")
legend(fig_leg)
%}

%{
%Total Spike Histograms
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = 0.512;
ex_c_idx = coherences == ex_c;
figure;
hold on
axs = [];
binEdges = 0:0.5:60;
for brain = brains
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j<=length(pulse_amps);
        if pulse
            output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
        else
            output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
        end
        load(strcat(output_stimpath, "/decisions.mat"), ...
            "totspikes_g1", "totspikes_g2", "totspikes_ns", "totspikes_int")
        axs(j) = subplot(3, 1, j);
        histogram(totspikes_g1(ex_c_idx, :), 'BinEdges', binEdges);
        if pulse
            title(sprintf("%0.0fnA Pulse Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
        else
            title(sprintf("%0.0fnA Galvanic Stimulation, %0.1f%% Coherence", [stim_amp*1e9, ex_c*100]))
        end
        if j == length(stim_amps)
            xlabel("Average Firing Rate (Hz)")
        elseif j == 2
            ylabel("Number of Neurons")
        end
    end
    linkaxes(axs);
    hold off
end
%}

function [time_idx, neuron_idx, g_idx] = get_spike_idx(g_spikes)
    [~, g_idx] = sort(sum(g_spikes, 1), 'descend');
    g_spikes = g_spikes(:, g_idx);
    [time_idx, neuron_idx] = find(g_spikes);
end
%}