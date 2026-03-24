%% --- MASTER SCRIPT: PARALLEL CALIBRATION GRID (SAFE VERSION) ---
clear; clc; close all;

tW_in = 15; tc_in = 50; 
dists_vec = [10, 20, 30, 40, 50, 60, 70,80];

B0_true_range = 0.1:0.1:0.9;      
Le_inv_true_range = 0.02:0.01:0.07; 
[B_mesh, L_mesh] = meshgrid(B0_true_range, Le_inv_true_range);
B_list = B_mesh(:); L_list = L_mesh(:);
num_sim = length(B_list);

results_map = zeros(num_sim, 4);

if isempty(gcp('nocreate')), parpool; end
fprintf('--- STARTING PARALLEL CALIBRATION (%d simulations) ---\n', num_sim);

parfor i = 1:num_sim

fprintf('Starting iteration %d of %d...\n', i, num_sim);

    % --- SETUP WORKER ---
    b_true = B_list(i);
    le_true = L_list(i);
    
    p_local = struct();
    p_local.Npart = 800000; p_local.dt = 0.2; p_local.tmax = 70;
    p_local.xs = 0; p_local.ys = 0; p_local.zs = 20; p_local.deltar = 0.7;
    p_local.A0 = 3.5; p_local.A1 = 3.5; p_local.A2 = 4.6;
    p_local.aa1 = 0; p_local.aa2 = 40; p_local.b1 = 0; p_local.b2 = 5;
    p_local.g2 = 0;
    p_local.g0 = b_true * le_true;
    p_local.etai = le_true - p_local.g0;
    p_local.g1 = p_local.g0;

    file_rand = 'random_math.dat'; 
    file_raw_worker = sprintf('scattering_worker_%d.dat', i);
    file_histo_worker = sprintf('histo_worker_%d.dat', i);

    % --- FASE 1: Simulazione Monte Carlo ---
    try
        scattering_layer(p_local, file_rand, file_raw_worker);
        
        % --- CONTROLLO ROBUSTO ESISTENZA FILE ---
        % Aspetta finché il file non è disponibile (max 5 secondi)
        timeout = 5; t_start = tic;
        while ~exist(file_raw_worker, 'file') && toc(t_start) < timeout
            pause(0.2); 
        end
        
        if exist(file_raw_worker, 'file')
            histogram_toro(p_local, dists_vec, file_raw_worker, file_histo_worker);
        else
            error('Timeout: file %s non creato dal worker %d', file_raw_worker, i);
        end
        
        % --- FASE 2: Estrazione Energie ---
        SimResults = extract_windows(p_local, dists_vec, file_histo_worker, tW_in, tc_in);
        
        % --- FASE 3: GRID SEARCH (Fit UHS) ---
        best_res = inf; B0_app = 0; Le_app = 0;
        p_fit = struct('v', p_local.A1, 'Fs', 1/p_local.dt, 'deltaT', tW_in, 't_coda', tc_in, 't_coda_in', tc_in); % Aggiunto t_coda_in se serve

        for b_test = 0.1:0.1:0.9 % Griglia leggermente più larga per velocità
            for le_test = 0.01:0.02:0.1
                Eth1 = computeEnergy3([b_test, le_test], dists_vec, 1, p_fit);
                Eth2 = computeEnergy3([b_test, le_test], dists_vec, 2, p_fit);
                Eth3 = computeEnergy3([b_test, le_test], dists_vec, 3, p_fit);
                
                res = sum((SimResults(:,2) - Eth1).^2) + ...
                      sum((SimResults(:,3) - Eth2).^2) + ...
                      sum((SimResults(:,4) - Eth3).^2);
                
                if res < best_res
                    best_res = res; B0_app = b_test; Le_app = le_test;
                end
            end
        end
        results_map(i, :) = [b_true, le_true, B0_app, Le_app];

    catch ME
        fprintf('Errore nel Worker %d: %s\n', i, ME.message);
    end
    
    % Pulizia finale
    % if exist(file_raw_worker, 'file'), delete(file_raw_worker); end
    % if exist(file_histo_worker, 'file'), delete(file_histo_worker); end
end

% CANCELLAZIONE MASSIVA ALLA FINE (Fuori dal parfor)
fprintf('Cleaning temporary files...\n');
delete('scattering_worker_*.dat');
delete('histo_worker_*.dat');

save('calibration_results.mat', 'results_map');
fprintf('--- CALIBRATION COMPLETED SUCCESSFULLY ---\n');

% save('calibration_results.mat', 'results_map');
% fprintf('--- CALIBRATION COMPLETED ---\n');

%%
% Dati della tua mappa (ho rimosso la riga degli zeri)
% 1. Caricamento dati dalla simulazione parallela
if exist('calibration_results.mat', 'file')
    load('calibration_results.mat'); % Carica la variabile 'results_map'
else
    error('File calibration_results.mat non trovato!');
end

% 2. Pulizia automatica: tiene solo le righe dove NON ci sono zeri
% (any(results_map, 2) restituisce vero se la riga ha almeno un elemento != 0)
clean_map = results_map(all(results_map ~= 0, 2), :);

fprintf('Simulazioni valide trovate: %d\n', size(clean_map, 1));

B0_t = clean_map(:,1); Le_t = clean_map(:,2); % Veri (HCTM)
B0_a = clean_map(:,3); Le_a = clean_map(:,4); % Apparenti (UHS)

% Calcolo degli scostamenti (Delta)
dB0 = B0_a - B0_t;
dLe = Le_a - Le_t;

% Costruzione matrice G (2 equazioni per ogni punto -> 16 righe totali)
N = length(B0_t);
G = zeros(2*N, 5);
d = zeros(2*N, 1);

for i = 1:N
    % Equazione per B0
    G(2*i-1, :) = [1, 2*B0_a(i), Le_a(i), 0, 0];
    d(2*i-1)    = dB0(i);
    % Equazione per Le^-1
    G(2*i, :)   = [0, 0, B0_a(i), 1, 2*Le_a(i)];
    d(2*i)      = dLe(i);
end

% Risoluzione minimi quadrati per trovare i coefficienti a1...a5
m_coeff = G \ d;

fprintf('Coefficienti calcolati:\n');
fprintf('a1: %.4f, a2: %.4f, a3: %.4f, a4: %.4f, a5: %.4f\n', m_coeff);
%%
% Esempio: i tuoi dati di Visso (inserisci i tuoi valori qui)
B0_visso_UHS = 0.69; 
Le_inv_visso_UHS = 0.073;

% Recupero coefficienti
a1=m_coeff(1); a2=m_coeff(2); a3=m_coeff(3); a4=m_coeff(4); a5=m_coeff(5);

% Calcolo del Bias (l'errore dovuto al mantello)
delta_B0 = a1 + 2*a2*B0_visso_UHS + a3*Le_inv_visso_UHS;
delta_Le = a4 + 2*a5*Le_inv_visso_UHS + a3*B0_visso_UHS;

% Risultato Finale: Parametri reali della crosta!
B0_corrected = B0_visso_UHS - delta_B0;
Le_inv_corrected = Le_inv_visso_UHS - delta_Le;

fprintf('\nParametri Corretti per Visso:\n');
fprintf('B0: %.3f, Le^-1: %.4f\n', B0_corrected, Le_inv_corrected);