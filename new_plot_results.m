clear;
sim_name = "EMBC";
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
default_colors = get(gca, "colororder");
start_trial = 1;
end_trial = 36;
num_trials = end_trial - start_trial + 1;
brains = 1:10;
num_brains = length(brains);
pulse_coherences = [-100, -78.8, -75.6, -72.4, -69.2, -66, -51.2, -25.6, 0, 25.6] / 100;
control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6] / 100;
galvanic_coherences = [-100, -51.2 -42.6, -39.4, -36.2, -33, -29.8, -25.6, 0, 25.6] / 100;
pulse_amps = [-10*1e-6];
dc_amps = [-28, 0]*1e-9;
stim_amps = [pulse_amps, dc_amps];
ex_c = 25.6/100;
ex_trial = 1;
ex_brain = 1;

%{
%Example FRs
figure;
axs = zeros(length(stim_amps), 1);
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    if pulse
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
            [sim_name, ex_brain, stim_amp*1e9]);
    elseif stim_amp == 0
        output_stimpath = sprintf("Simulation %s/brain1/data/%0.1fnA_galvanic", ...
            [sim_name, stim_amp*1e9]);
    else
        output_stimpath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
            [sim_name, ex_brain, stim_amp*1e9]);
    end
    try
        load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), "pop_frs")
    catch
        continue
    end
    axs(j) = subplot(length(stim_amps), 1, j);
    hold on
    for i = 1:p+2
        plot(t, pop_frs(:, i))
    end
    xline(t_task, "g--")
    xline(t_taskoff, "r--")
    hold off
    if pulse
        title("Pulsatile Stimulation")
    elseif stim_amp == 0
        title("Control")
    else
        title("Galvanic Stimualtion")
    end
    if j == length(stim_amps)
        xlabel("Time (s)")
    elseif j == ceil(length(stim_amps)/2)
        ylabel("Population Firing Rate (Hz)")
    elseif j == 1
        legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory", "Stimulation On", "Stimulation Off"], ...
            "Location", "northwest", 'FontSize', 6)
    end
end
linkaxes(axs)
%}

%{
%Example Raster: P1 vs. Distance
load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, ex_brain]), "ball_r")
p1_only = figure;
top_N = 12;
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
    spikes = zeros(length(t), top_N);
    neuron_num = 1:N;
    for nn = 1:top_N
        for spike_idx = recspikes(int2str(nn))
            spikes(spike_idx, nn) = 1;
        end
    end
    
    [g1_time_idx, g1_neuron_idx, g1_idx] = get_spike_idx(spikes);
    
    figure(p1_only);
    subplot(length(stim_amps), 1, j)
    scatter(t(g1_time_idx), ball_r(g1_idx(g1_neuron_idx))*1e6, "Marker", "|", ...
        "MarkerFaceColor", default_colors(1, :), "MarkerEdgeColor", default_colors(1, :))
    if j == length(stim_amps)
        xlabel("Time (s)")
    elseif j == 2
        ylabel("Distance from Electrode (um)")
    end
    if pulse
        title("Pulsatile Stimulation")
    else
        if abs(stim_amp) > 0
            title("Galvanic Stimulation")
        else
            title("Control")
        end
    end
end
%}

%{
%Example ISI
load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, ex_brain]), "ball_r")
p1_vs_distance = figure;
p1_total = figure;
ax = [];
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
    var_isi = zeros(num_group, 1);
    isi = [];
    for nn = 1:num_group
        spikeidx = recspikes(int2str(nn));
        var_isi(nn) = var(t(diff(spikeidx)));
        isi = [isi, t(diff(spikeidx))];
    end
    
    figure(p1_vs_distance);
    ax(j) = subplot(length(stim_amps), 1, j);
    scatter(ball_r, var_isi)
    if pulse
        title("Pulsatile Stimulation")
    else
        if abs(stim_amp) > 0
            title("Galvanic Stimulation")
        else
            title("Control")
        end
    end
    var(isi)
    figure(p1_total)
    hold on
    bar(var(isi))
end
hold off
figure(p1_vs_distance)
linkaxes(ax)
%}

%{
%Single Neuron FRs
load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, ex_brain]), "ball_r")
g1_taskfrs = zeros(num_group, 3);
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
    for nn = 1:num_group
        spiketimes = t(recspikes(int2str(nn)));
        g1_taskfrs(nn, j) = sum(spiketimes>=(t_task+1.5) & spiketimes<t_taskoff) * 2; %t_task + 1.5 to minimize ramping effect
    end
end
top_N = num_group;
figure;
hold on
scatter(ball_r(1:top_N)*1e6, g1_taskfrs(1:top_N, 1), [], ones(top_N, 3).*default_colors(7, :), 'filled')
scatter(ball_r(1:top_N)*1e6, g1_taskfrs(1:top_N, 2), [], ones(top_N, 3).*default_colors(5, :), 'filled')
scatter(ball_r(1:top_N)*1e6, g1_taskfrs(1:top_N, 3), [], "k", 'filled')
hold off
xlabel("Distance from Electrode (um)")
ylabel("End of Task Firing Rate (Hz)")
legend(["Pulsatile", "Galvanic", "Control"])
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


%Final Accuracies
pulse_acc = zeros(length(brains), length(pulse_coherences));
galvanic_acc = zeros(length(brains), length(galvanic_coherences));
ctrl_acc = zeros(1, length(control_coherences));
c = -1:0.01:1;
for brain = brains
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
            stim_coherences = pulse_coherences;
        elseif stim_amp == 0 %Control
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
            stim_coherences = control_coherences;
            if brain ~= 1
                continue
            end
        else
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
            stim_coherences = galvanic_coherences;
        end
        load(strcat(datapath, "/decisions.mat"), "avg_final_acc", "final_decisions");
        [nodec_trial, nodec_c] = find(~final_decisions);
        for nodec = 1:length(nodec_trial)
            fprintf("Brain %0.0f, c=%0.3f, trial = %0.0f \n", ...
                [brain, stim_coherences(nodec_c(nodec)), nodec_trial(nodec)])
        end
        if pulse
            pulse_acc(brain, :) = avg_final_acc;
        elseif stim_amp == 0
            ctrl_acc = avg_final_acc;
        else
            galvanic_acc(brain, :) = avg_final_acc;
        end
    end
end

figure;
hold on
errorbar(pulse_coherences, mean(pulse_acc, 1), std(pulse_acc, [], 1)/sqrt(10), 'Color', default_colors(7, :))
errorbar(galvanic_coherences, mean(galvanic_acc, 1), std(galvanic_acc, [], 1)/sqrt(10), 'Color', default_colors(5, :))
plot(control_coherences, ctrl_acc, 'ko-')
scatter(pulse_coherences, pulse_acc', [], default_colors(7, :).*ones(length(pulse_acc), 3))
scatter(galvanic_coherences, galvanic_acc', [], default_colors(5, :).*ones(length(pulse_acc), 3))
hold off
xlabel("Coherence (%)")
ylabel("% of trials P1 wins")
legend("Pulsatile", "Galvanic", "Control")
%}


%Final Decision Times
pulse_dt = zeros(length(brains), length(pulse_coherences));
galvanic_dt = zeros(length(brains), length(galvanic_coherences));
ctrl_dt = zeros(1, length(control_coherences));
c = -1:0.01:1;
for brain = brains
    for j = 1:length(stim_amps)
        stim_amp = stim_amps(j);
        pulse = j <= length(pulse_amps);
        if pulse
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_pulse", ...
                [sim_name, brain, stim_amp*1e9]);
        elseif stim_amp == 0 %Control
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
            if brain ~= 1
                continue
            end
        else
            datapath = sprintf("Simulation %s/brain%0.0f/data/%0.1fnA_galvanic", ...
                [sim_name, brain, stim_amp*1e9]);
        end
        load(strcat(datapath, "/decisions.mat"), "final_decision_times", "final_decisions");
        size(final_decisions)
        if pulse
            pulse_dt(brain, :) = mean(final_decision_times, 1, 'omitnan');
        elseif stim_amp == 0
            ctrl_dt(brain, :) = mean(final_decision_times, 1, 'omitnan');
        else
            galvanic_dt(brain, :) = mean(final_decision_times, 1, 'omitnan');
        end
    end
end

figure;
hold on
errorbar(pulse_coherences, mean(pulse_dt, 1), std(pulse_dt, [], 1)/sqrt(num_brains), 'Color', default_colors(7, :))
errorbar(galvanic_coherences, mean(galvanic_dt, 1), std(galvanic_dt, [], 1)/sqrt(num_brains), 'Color', default_colors(5, :))
plot(control_coherences, ctrl_dt, 'ko-')
scatter(pulse_coherences, pulse_dt', [], default_colors(7, :).*ones(length(pulse_dt), 3))
scatter(galvanic_coherences, galvanic_dt', [], default_colors(5, :).*ones(length(pulse_dt), 3))
hold off
xlabel("Coherence (%)")
ylabel("Decision Time (s)")
legend("Pulsatile", "Galvanic", "Control")
%}

%{
%Spike Timing Correlation 
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = 0.256;
ex_brain = 1;
timeidx = t >= t_task & t < t_taskoff;
load(sprintf("Simulation %s/brain%0.0f/r.mat", [sim_name, ex_brain]), "ball_r")
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
    fprintf("Stimulation %0.1fnA \n", stim_amp*1e9)
    neuron_frs = spikes2neuron_frs(spikes, dt);
    
    g1_frs = neuron_frs(timeidx, 1:num_group);
    g1_frs = g1_frs(:, ~all(g1_frs==0, 1)); %omit neurons that never fire
    g1_corr = cov(g1_frs) ./ sqrt(var(g1_frs) .* var(g1_frs)');
    g1_corr = g1_corr(~eye(size(g1_corr))); %omit diagonal
    g1_avg_corr = mean(g1_corr);
    fprintf("g1 : %0.4f \n", g1_avg_corr)
    
    g2_frs = neuron_frs(timeidx, num_group+1:2*num_group);
    g2_frs = g2_frs(:, ~all(g2_frs==0, 1)); %omit neurons that never fire
    g2_corr = cov(g2_frs) ./ sqrt(var(g2_frs) .* var(g2_frs)');
    g2_corr = g2_corr(~eye(size(g2_corr))); %omit diagonal
    g2_avg_corr = mean(g2_corr);
    fprintf("g2 : %0.4f \n", g2_avg_corr)
    
    ns_frs = neuron_frs(timeidx, 2*num_group+1:N_E);
    ns_frs = ns_frs(:, ~all(ns_frs==0, 1)); %omit neurons that never fire
    ns_corr = cov(ns_frs) ./ sqrt(var(ns_frs) .* var(ns_frs)');
    ns_corr = ns_corr(~eye(size(ns_corr))); %omit diagonal
    ns_avg_corr = mean(ns_corr);
    fprintf("ns : %0.4f \n", ns_avg_corr)
    
    int_frs = neuron_frs(timeidx, N_E+1:end);
    int_frs = int_frs(:, ~all(int_frs==0, 1)); %omit neurons that never fire
    int_corr = cov(int_frs) ./ sqrt(var(int_frs) .* var(int_frs)');
    int_corr = int_corr(~eye(size(int_corr))); %omit diagonal
    int_avg_corr = mean(int_corr);
    fprintf("int : %0.4f \n", int_avg_corr)
    
end
%}

function [time_idx, neuron_idx, g_idx] = get_spike_idx(g_spikes)
    [~, g_idx] = sort(sum(g_spikes, 1), 'descend');
    g_spikes = g_spikes(:, g_idx);
    [time_idx, neuron_idx] = find(g_spikes);
end

function neuron_frs = spikes2neuron_frs(spikes, dt)
    win_size = 5e-3;
    avg_win_size = 50e-3;
    w = ones(floor(win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, spikes) ./ dt;
    w = ones(floor(avg_win_size/dt), 1);
    w = w ./ length(w);
    neuron_frs = filter(w, 1, neuron_frs);
end
%}