clear;
sim_num = 2;
sim_path = sprintf("Simulation %0.0f", sim_num);
load(strcat(sim_path, "/bam_constants.mat"))
dc_colors = ["g", "k", "r"];

%{
%Example FRs
ex_c = -0.032;
ex_trial = 1;
%ex_stim_amp = -1*1e-9;
pulse = false;
%ex_stim_amp = -100*1e-9;
ex_stim_amp = 0;
figure;
hold on
if pulse
    output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", [sim_num, ex_stim_amp*1e9]);
    title(sprintf("%0.0fnA Pulse Stimulation, Coherence: %0.1f", [ex_stim_amp*1e9, ex_c*100]))
else
    output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
        [sim_num, ex_stim_amp*1e9]);
    title(sprintf("%0.0fnA Galvanic Stimulation, Coherence: %0.1f", [ex_stim_amp*1e9, ex_c*100]))
end
load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), "pop_frs")
for i = 1:p+2
    plot(t, pop_frs(:, i))
end
hold off
legend([compose("Selective Population %0.0f", 1:p), "Non-Selective", "Inhibitory"])
%}

%{
%Example FRs with Raster
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = -0.032;
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    figure;
    hold on
    if pulse
        output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", [sim_num, stim_amp*1e9]);
        title(sprintf("%0.0fnA Pulse Stimulation", stim_amp*1e9))
    else
        output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
            [sim_num, stim_amp*1e9]);
        title(sprintf("%0.0fnA Galvanic Stimulation", stim_amp*1e9))
    end
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
        "pop_frs", "recspikes")
    axon_SA = 100*(1e-6^2);
    ball_r = sqrt(axon_SA ./ (abs(electric_r)*4*pi));
    neuron_num = 1:N;
    for nn = 1:num_group
        spike_times = t(recspikes(int2str(nn)));
        scatter(spike_times, ones(size(spike_times))*ball_r(nn)*1e6, 'Marker', '|')
    end
    xlabel("Time (s)")
    ylabel("Distance From Electrode (um)")
end
hold off
%}

%{
%Example Raster by pop-type and total activity 
ex_trial = 1;
default_colors = get(gca, 'colororder');
ex_c = -0.032;
for j = 1:length(stim_amps)
    stim_amp = stim_amps(j);
    pulse = j<=length(pulse_amps);
    figure;
    hold on
    if pulse
        output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_pulse", [sim_num, stim_amp*1e9]);
        title(sprintf("%0.0fnA Pulse Stimulation", stim_amp*1e9))
    else
        output_stimpath = sprintf("Simulation %0.0f/data/%0.0fnA_galvanic", ...
            [sim_num, stim_amp*1e9]);
        title(sprintf("%0.0fnA Galvanic Stimulation", stim_amp*1e9))
    end
    load(strcat(output_stimpath, sprintf("/c=%0.3f/trial%0.0f.mat", [ex_c, ex_trial])), ...
        "pop_frs", "recspikes")
    axon_SA = 100*(1e-6^2);
    ball_r = sqrt(axon_SA ./ (abs(electric_r)*4*pi));
    neuron_num = 1:N;
    for nn = 1:N
        if nn <= num_group
            color_idx = 1;
        elseif nn > num_group && nn <= 2*num_group
            color_idx = 2;
        elseif nn > 2*num_group && nn < N_E
            color_idx = 3;
        else
            color_idx = 4;
        end
        spike_times = t(recspikes(int2str(nn)));
        scatter(spike_times, ones(size(spike_times))*length(spike_times), 'Marker', '|', 'MarkerEdgeColor', default_colors(color_idx, :))
    end
    xlabel("Time (s)")
    ylabel("Total Number of spikes")
end
hold off
%}
