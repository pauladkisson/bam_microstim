%clear;
<<<<<<< HEAD
sim_name = "%Activation Equivalence Disconnected";
=======
sim_name = "Test";
>>>>>>> 50d200f4fda53668a87df9c65f3f4f6722acd698
sim_path = sprintf("Simulation %s", sim_name);
load(strcat(sim_path, "/bam_constants.mat"))
figure;
default_colors = get(gca, "colororder");
start_trial = 1;
end_trial = 36;
num_trials = end_trial - start_trial + 1;
%brains = 1:10;
brains = 1:3;
num_brains = length(brains);
num_batch = 3;

pulse_coherences = [-100, -78.8, -75.6, -72.4, -69.2, -66, -51.2, -25.6, 0, 25.6] / 100;
control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6] / 100;
%galvanic_coherences = [-100, -51.2 -42.6, -39.4, -36.2, -33, -29.8, -25.6, 0, 25.6] / 100;
galvanic_coherences = [-100, -75.6, -51.2, -25.6, 0, 25.6] / 100;
%}
%{
pulse_coherences = [0];
control_coherences = [0];
galvanic_coherences = [0];
%}
pulse_amps = [-10*1e-6];
<<<<<<< HEAD
dc_amps = [-120, 0]*1e-9;
=======
dc_amps = [-110, 0]*1e-9;
>>>>>>> 50d200f4fda53668a87df9c65f3f4f6722acd698
stim_amps = [pulse_amps, dc_amps];

%{
ex_c = 0/100;
ex_trial = 1;
ex_brain = 1;
ex_stim_j = 1;
plot_name = "single_stim"; % or 'subplot' or 'p1_only'
plot_frs(sim_name, pulse_amps, stim_amps, p, t, default_colors, ex_stim_j, ex_brain, ex_c, ex_trial, plot_name)
%}

%{
top_N = num_group*0.1;
ex_neurons = [5, 6];
plot_name = "grouped_stim"; % or 'single_stim' or 'grouped_stim'
plot_rasters(sim_name, pulse_amps, stim_amps, ex_neurons, t, t_task, t_taskoff, stim_freq, default_colors, top_N, ex_stim_j, ex_brain, ex_c, ex_trial, plot_name)
%}

%{
win_size = floor(0.250 / dt); %250ms moving window
%cv_window = t >= 2.5 & t<3; %Plotting window
cv_window = t >= t_task & t<t_taskoff; %Plotting window
ex_neuron = 7;
top_N = floor(num_group);
plot_name = "p1_wins"; % or 'ex_trial' or 'p1_wins'
plot_cv(sim_name, pulse_amps, stim_amps, t, N, top_N, num_group, ...
                 win_size, cv_window, default_colors, ex_brain, ex_c, ex_trial, ...
                 ex_neuron, brains, num_brains, pulse_coherences, galvanic_coherences, control_coherences, ...
                 start_trial, end_trial, num_trials, plot_name)
%}

%{
idx_diff = stim_ind+1;% how far off timing is from pulse timing + 1 to account for t(1) = 0
plot_phaselock(pulse_amps, stim_amps, t, t_task, t_taskoff, stim_freq, num_group, ...
                        idx_diff, default_colors, brains, num_brains, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials)
%}

%{
N_start = 1;
N_end = floor(num_group);
win_start = 2.5;
win_stop = 3;
c_win = 300*1e-6;
c = 0;
brains = [1, 2];
num_brains = 2;
sim_names = ["EMBC I_b100", "EMBC Disconnected"];
plot_sync(sim_names, pulse_amps, stim_amps, t, num_group, ...
                        brains, num_brains, N_start, N_end, ...
                        win_start, win_stop, c_win, c, ...
                        pulse_coherences, galvanic_coherences, control_coherences, ...
                        start_trial, end_trial, num_trials, default_colors)
%}


win_start = t_task + stim_ind*dt; % to account for onset spike of pulse
win_stop = t_task + 0.1;
ex_c = 0;
%  plot_name = 'ex_c' or 'p1_wins' or 'p1_loses'
plot_name = "p1_wins";
plot_frdist(sim_name, ex_c, pulse_amps, stim_amps, t, num_group, win_start, ...
                     win_stop, default_colors, brains, num_brains, ...
                     pulse_coherences, galvanic_coherences, control_coherences, ...
                     start_trial, end_trial, num_trials, plot_name);
%}

%{
plot_decisions(sim_name, pulse_amps, stim_amps, default_colors, brains, ...
               num_brains, num_batch, ...
               pulse_coherences, galvanic_coherences, control_coherences)
%}
