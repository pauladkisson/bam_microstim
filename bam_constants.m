%%% Paul Adkisson 
%%% 9.2.21
%%% Purpose Define Constants for Wang 2002 Biophysical Attractor
%%% Model (BAM)
%%% Adds pulse-pulse, pulse-spontaneous, spontaneous-pulse blocking
clear;
tic;

%% Simulation Parameters
sim_name = "Test";
sim_path = sprintf("Simulation %s", sim_name);
mkdir(sim_path)
dt = 0.05e-3; %ms
t_span = 4;
t = 0:dt:t_span;
start_trial = 1;
end_trial = 12;

%% Network Parameters
%percent_size = 0.5;
percent_size = 0.15;
N_E = floor(1600 * percent_size);
%N_I = floor(400 * percent_size);
N_I = 0;
N = N_E + N_I;
%f = 0.15; % Fraction of cortical neurons asctivated by one type of stimulus
f = 0.5;
p = 2; % Number of different types of stimuli
num_selective = floor(p*f*N_E);
num_group = floor(f*N_E);
%w_plus = 1.7; % Strength of "strong" synapses in the BAM network
%w_minus = 1 - f*(w_plus - 1)/(1-f); %Strength of "weak" synapses in BAM
%w = 1; %Strength of normal synapses in BAM
w_plus = 0;
w_minus = 0;
w = 0;
start_brain = 1;
end_brain = 1;
brains = start_brain:end_brain;
GenerateBAM(brains, N_E, N_I, f, p, w_plus, w_minus, w, sim_path);
GenerateConductances(N_E, N_I, sim_path)
pop_type = ones(N, 1);
pop_type(N_E+1:end) = 2; % population_type = 1 for pyr, 2 for int

%% Input Parameters
fr_bg = 2400;
% Synaptic Conductance = [pyramidal, interneuron]
G_ampa_ext = [2.1, 1.62]*1e-9; %nS
%pulse_coherences = [-100, -78.8, -75.6, -72.4, -69.2, -66, -51.2, -25.6, 0, 25.6] / 100;
%control_coherences = [-100, -51.2, -25.6, -12.8, -6.4, -3.2, 0, 3.2, 6.4, 12.8, 25.6] / 100;
%galvanic_coherences = [-100, -51.2 -42.6, -39.4, -36.2, -33, -29.8, -25.6, 0, 25.6] / 100;
pulse_coherences = [0] / 100;
control_coherences = [0] / 100;
galvanic_coherences = [0] / 100;
coherences = union(union(pulse_coherences, galvanic_coherences), control_coherences, 'sorted');
max_fr_task = 80;
t_task = 1;
t_taskoff = 3;
GenerateSpikes(fr_bg, max_fr_task, coherences, f, N_E, N_I, t_task, ...
    t_taskoff, t, start_trial, end_trial, sim_path);

%% LIF Parameters
% Parameter = [pyramidal, interneuron]
C = [0.5, 0.2]*1e-9; %nF
gL = [25, 20]*1e-9; %nS
EL = -70e-3; %mV
Vs = -50e-3; %mV
Vr = -55e-3;
tau_r = [2, 1]*1e-3; %ms
syn_delay = 0.5e-3; %ms
refract_ind = floor(tau_r/dt);
delay_ind = floor(syn_delay/dt);
tau_AMPA = 2e-3; %ms
tau_NMDA_1 = 2e-3; %ms
tau_NMDA_2 = 100e-3; %ms
tau_GABA = 5e-3; %ms
alpha = 500; %Hz

%% Microstimulation Parameters
stim_duration = 300e-6; %us / phase
stim_ind = floor(stim_duration*2 / dt);
stim_freq = 200; %Hz
depol_block_thresh = 800*1e-12;
depol_block_factor = 3;
pulse_amps = [-5*1e-6];
dc_amps = [-40, 0]*1e-9;
stim_amps = [pulse_amps, dc_amps];
GenerateMicroStim(t, t_task, t_taskoff, stim_duration, stim_freq, ...
                  depol_block_thresh, depol_block_factor, pulse_amps, ...
                  dc_amps, N, num_group, brains, sim_path);

%% Firing Rate Parameters
win_size = 5e-3;
win_index = floor(win_size / dt);
avg_win_size = 50e-3;
avg_win_index = floor(avg_win_size / dt);

%% t_b
t_b = [0 0 0 
    4.5000         0   25.0000

    9.0000         0   25.0000

    13.5000         0   25.0000

    18.0000         0   25.0000

    22.5000         0   25.0000

    27.0000         0   25.0000

    31.5000         0   20.0000

    36.0000         0   20.0000

    40.5000         0   20.0000

    45.0000         0   20.0000

    49.5000         0   14.5000

    54.0000         0   13.40

    58.5000         0   12.6

    63.0000         0    10.250

    67.5000         0    8.6000

    72.0000         0    7.0000

    76.5000         0    6.1000

    81.0000         0    5.8000

    85.5000         0    5.5000

    90.0000         0    5.4000

    96.0000         0    4.8000

    108.0000         0    4.7600

    120.0000         0    4.7500

    132.0000         0    4.7100

    144.0000         0    4.7100

    156.0000         0    4.7300

    168.0000         0    4.7300

    180.0000         0    4.7300

    192.0000         0    4.7300

    204.0000         0    4.7300

    216.0000         0    5.1502

    228.0000         0    5.1852

    240.0000         0    5.5000

    252.0000         0    6.2000

    264.0000         0    7.0000

    276.0000         0    7.7000

    288.0000         0    8.3000

    300.0000         0   16.2580

    312.0000         0   31.6985

    324.0000         0   68.5536

    336.0000         0  156.5234

    348.0000         0  366.4995

    360.0000         0  867.6938
    
    500.0000         0  867.6938];

load('params_tPP_tPS_tSP_from_10_27_20.mat')
best_tref_per_c = best_tref_per_c(1:size(t_b, 1), :);
t_pp = best_tref_per_c(:, 1)*1e-3;
t_ps = best_tref_per_c(:, 2)*1e-3;
[min_ps, min_ps_idx] = min(t_ps(2:end));
t_ps(1:min_ps_idx) = 0; %LIF model already accounts for low-amplitude insufficiency
t_sp = best_tref_per_c(:, 3)*1e-3;
I_b = t_b(:, 1)*100*1e-9/360;
 
%% Save
save_path = strcat(sim_path, "/bam_constants.mat");
save(save_path)
toc


