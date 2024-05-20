%% Clear the environment and the command line
clear;
clc;
close all;

%% Define the input parmeters
% ---------------------------------------------------------
num_servers = 2;        % cannot be changed in this code
% ---------------------------------------------------------

mean_arr_rate = 2;        % Mean arrival rate (cust/min)
mean_ser_rate = 1.5;      % Mean service rate (cust/min)
num_planes_landing = 5000; % Number of landing planes
num_planes_takeoff = 5000; % Number of takeoff planes     

% Variable initializations
%--landing--
service_start_time_landing = zeros(1, num_planes_landing);
service_time_landing = zeros(1, num_planes_landing);
service_end_time_landing = zeros(1, num_planes_landing);
%--takeoff--
service_start_time_takeoff = zeros(1, num_planes_takeoff);
service_time_takeoff= zeros(1, num_planes_takeoff);
service_end_time_takeoff= zeros(1, num_planes_takeoff);

%% Calculation Module
assert(num_servers == 2);

% Define random arrival time distribution

%--landing--
pd_arr = makedist('exponential', 'mu', 1/mean_arr_rate); %Exponential Distribution Object Created Using makedist('Exponential')
rand_arr_time_landing = random(pd_arr, 1, num_planes_landing);
cum_arr_time_landing = cumsum(rand_arr_time_landing);   % Arrival time in real scale

%--takeoff--
rand_arr_time_takeoff= random(pd_arr, 1,num_planes_takeoff);
cum_arr_time_takeoff= cumsum(rand_arr_time_takeoff);   % Arrival time in real scale

% Define random service time distribution when only one server is utilized
pd_first_ser = makedist('exponential', 'mu', 1/mean_ser_rate);
rand_f_ser_time_landing  = random(pd_first_ser, 1, num_planes_landing); 
rand_f_ser_time_takeoff= random(pd_first_ser, 1, num_planes_takeoff); 

% Define random service time distribution when both servers are utilized
pd_ser = makedist('exponential', 'mu', 1/(num_servers * mean_ser_rate));
rand_ser_time_landing = random(pd_ser, 1, num_planes_landing);
rand_ser_time_takeoff= random(pd_ser, 1, num_planes_takeoff);

%% Simulation

% Setup the first landing plane.
land = 1;
service_start_time_landing(land) = cum_arr_time_landing(land);
service_time_landing(land) = rand_f_ser_time_landing (land);
service_end_time_landing(land) = service_start_time_landing(land) + service_time_landing(land);

% Setup the first takeoff plane.
takeoff = 1;
service_start_time_takeoff(takeoff) = cum_arr_time_takeoff(takeoff);
service_time_takeoff(takeoff) = rand_f_ser_time_takeoff(takeoff);
service_end_time_takeoff(takeoff) = service_start_time_takeoff(takeoff) + service_time_takeoff(takeoff);

% Setup the second landing plane
land = 2;
service_start_time_landing(land) = cum_arr_time_landing(land);
% New plane arrives when one server is occupied
if cum_arr_time_landing(land) < service_end_time_landing(land - 1)
    service_time_landing(land) = rand_ser_time_landing(land);
% New plane arrives when both servers are free
else
    service_time_landing(land) = rand_f_ser_time_landing (land);
end
service_end_time_landing(land) = service_start_time_landing(land) + service_time_landing(land);

% Setup the second takeoff plane
takeoff = 2;
service_start_time_takeoff(takeoff) = cum_arr_time_takeoff(takeoff);
% New plane arrives when one server is occupied
if cum_arr_time_takeoff(takeoff) < service_end_time_takeoff(takeoff- 1)
    service_time_takeoff(takeoff) = rand_ser_time_takeoff(takeoff);
% New plane arrives when both servers are free
else
    service_time_takeoff(takeoff) = rand_f_ser_time_takeoff(takeoff);
end
service_end_time_takeoff(takeoff) = service_start_time_takeoff(takeoff) + service_time_takeoff(takeoff);

% Calculation of Service Time for Remaining Planes

for land = 3:1:num_planes_landing
    % New plane arrives when both servers are occupied
    if cum_arr_time_landing(land) < service_end_time_landing(land - 2)
        service_start_time_landing(land) = service_end_time_landing(land - 2);
        service_time_landing(land) = rand_ser_time_landing(land);
    % New plane arrives when one server is occupied
    elseif cum_arr_time_landing(land) < service_end_time_landing(land - 1)
        service_start_time_landing(land) = service_end_time_landing(land - 1);
        service_time_landing(land) = rand_ser_time_landing(land);
    % New plane arrives when both servers are free
    else
        service_start_time_landing(land) = cum_arr_time_landing(land);
        service_time_landing(land) = rand_f_ser_time_landing (land);
    end
    service_end_time_landing(land) = service_start_time_landing(land) + service_time_landing(land);
end

for takeoff = 3:1:num_planes_takeoff
    % New plane arrives when both servers are occupied
    if cum_arr_time_takeoff(takeoff) < service_end_time_takeoff(takeoff- 2)
        service_start_time_takeoff(takeoff) = service_end_time_takeoff(takeoff- 2);
        service_time_takeoff(takeoff) = rand_ser_time_takeoff(takeoff);
    % New plane arrives when one server is occupied
    elseif cum_arr_time_takeoff(takeoff) < service_end_time_takeoff(takeoff - 1)
        service_start_time_takeoff(takeoff) = service_end_time_takeoff(takeoff - 1);
        service_time_takeoff(takeoff) = rand_ser_time_takeoff(takeoff);
    % New plane arrives when both servers are free
    else
        service_start_time_takeoff(takeoff) = cum_arr_time_takeoff(takeoff);
        service_time_takeoff(takeoff) = rand_f_ser_time_takeoff(takeoff);
    end
    service_end_time_takeoff(takeoff) = service_start_time_takeoff(takeoff) + service_time_takeoff(takeoff);
end
%% Post Simulation Computations

system_time_landing = service_end_time_landing - cum_arr_time_landing;
q_wait_time_landing = service_start_time_landing - cum_arr_time_landing;

system_time_takeoff = service_end_time_takeoff- cum_arr_time_takeoff;
q_wait_time_takeoff = service_start_time_takeoff- cum_arr_time_takeoff;

cum_q_wait_time_landing = cumsum(q_wait_time_landing);
run_mean_q_wait_time_landing = cum_q_wait_time_landing ./ (1:num_planes_landing);

cum_q_wait_time_takeoff = cumsum(q_wait_time_takeoff);
run_mean_q_wait_time_takeoff = cum_q_wait_time_takeoff ./ (1:num_planes_takeoff);

%% Calculate Queue Length at differnt points in time
delta_t = 1;                            % Time resolution is 1 minute

time_line_landing = 0:delta_t:floor(service_end_time_landing(end));
queue_len_landing = zeros(size(time_line_landing));

time_line_takeoff = 0:delta_t:floor(service_end_time_takeoff(end));
queue_len_takeoff = zeros(size(time_line_takeoff));

for time_ind = 1:size(time_line_landing, 2)
    for land_ind = 1:num_planes_landing
        if (cum_arr_time_landing(land_ind) < time_line_landing(time_ind) && ...
                service_end_time_landing(land_ind) > time_line_landing(time_ind) && ...
                service_start_time_landing(land_ind) > time_line_landing(time_ind))
            queue_len_landing(time_ind) = queue_len_landing(time_ind) + 1;
        end
    end
end
for time_ind = 1:size(time_line_takeoff, 2)
    for takeoff_ind = 1:num_planes_takeoff
        if (cum_arr_time_takeoff(takeoff_ind) < time_line_takeoff(time_ind) && ...
                service_end_time_takeoff(takeoff_ind) > time_line_takeoff(time_ind) && ...
                service_start_time_takeoff(takeoff_ind) > time_line_takeoff(time_ind))
            queue_len_takeoff(time_ind) = queue_len_takeoff(time_ind) + 1;
        end
    end
end
cum_queue_len_landing = cumsum(queue_len_landing);
run_mean_queue_len_landing = cum_queue_len_landing ./ (1:size(queue_len_landing, 2));

cum_queue_len_takeoff = cumsum(queue_len_takeoff);
run_mean_queue_len_takeoff = cum_queue_len_takeoff ./ (1:size(queue_len_takeoff, 2));
%% Analytic computations
% util_rate = mean_arr_rate / (num_servers * mean_ser_rate);
mean_q_wait_time_landing = run_mean_q_wait_time_landing(end);
mean_q_length_landing = run_mean_queue_len_landing(end);
mean_system_time_landing = mean(system_time_landing);

mean_q_wait_time_takeoff= run_mean_q_wait_time_takeoff(end);
mean_q_length_takeoff= run_mean_queue_len_takeoff(end);
mean_system_time_takeoff= mean(system_time_takeoff);

std_q_wait_time_landing = std(q_wait_time_landing);
std_q_length_landing = std(queue_len_landing);
std_system_time_landing = std(system_time_landing);

std_q_wait_time_takeoff = std(q_wait_time_takeoff);
std_q_length_takeoff = std(queue_len_takeoff);
std_system_time_takeoff = std(system_time_takeoff);

ste_q_wait_time_landing = std_q_wait_time_landing/sqrt(size(q_wait_time_landing, 2));
ste_q_length_landing = std_q_length_landing/sqrt(size(queue_len_landing, 2));
ste_system_time_landing = std_system_time_landing/sqrt(size(system_time_landing, 2));

ste_q_wait_time_takeoff = std_q_wait_time_takeoff/sqrt(size(q_wait_time_takeoff, 2));
ste_q_length_takeoff = std_q_length_takeoff/sqrt(size(queue_len_takeoff, 2));
ste_system_time_takeoff = std_system_time_takeoff/sqrt(size(system_time_takeoff, 2));

%% Print the computation results

fprintf('Mean Queue Waiting Time(Wqmc) for Landing Planes: %.2f minutes\n', mean_q_wait_time_landing)
fprintf('Mean Queue Length (Lqmc) for Landing Planes: %.2f planes \n', mean_q_length_landing)
fprintf('Mean System Time (Wmc) for Landing Planes: %.2f minutes\n\n', mean_system_time_landing)

fprintf('STD Queue Waiting Time (Wqmc)for Landing Planes: %.2f minutes\n', std_q_wait_time_landing)
fprintf('STD Queue Length (Lqmc) for Landing Planes: %.2f planes \n', std_q_length_landing)
fprintf('STD System Time (Wmc) for Landing Planes: %.2f minutes\n\n', std_system_time_landing)

fprintf('STE Queue Waiting Time (Wqmc) for Landing Planes: %.2f minutes\n', ste_q_wait_time_landing)
fprintf('STE Queue Length (Lqmc) for Landing Planes: %.2f planes \n', ste_q_length_landing)
fprintf('STE System Time (Wmc) for Landing Planes: %.2f minutes\n\n', ste_system_time_landing)


fprintf('Mean Queue Waiting Time (Wqmc) for Take-Off Planes: %.2f minutes\n', mean_q_wait_time_takeoff)
fprintf('Mean Queue Length (Lqmc) for Take-Off Planes: %.2f planes \n', mean_q_length_takeoff)
fprintf('Mean System Time (Wmc) for Take-Off Planes: %.2f minutes\n\n', mean_system_time_takeoff)

fprintf('STD Queue Waiting Time (Wqmc) for Take-Off Planes: %.2f minutes\n', std_q_wait_time_takeoff)
fprintf('STD Queue Length (Lqmc) for Take-Off Planes: %.2f planes \n', std_q_length_takeoff)
fprintf('STD System Time (Wmc) for Take-Off Planes: %.2f minutes\n\n', std_system_time_takeoff)

fprintf('STE Queue Waiting Time (Wqmc) for Take-Off Planes: %.2f minutes\n', ste_q_wait_time_takeoff)
fprintf('STE Queue Length (Lqmc) for Take-Off Planes: %.2f landing planes \n', ste_q_length_takeoff)
fprintf('STE System Time (Wmc) for Take-Off Planes: %.2f minutes\n\n', ste_system_time_takeoff)

%% Plot Outcomes

analytic_ss_queue_wait_time = 0.533;
analytic_ss_queue_len = 1.066;
%landing
figure(1)
plot(cum_arr_time_landing, run_mean_q_wait_time_landing)
hold on;
yline(analytic_ss_queue_wait_time, '--r');
hold off;
legend('MC Simulation Wq', 'Analytic SS Wq');
xlabel('Arrival time for Landing Planes')
ylabel('Mean Time in Queue for Landing Planes')
title('Mean Time in Queue Upon Arrival for Landing Planes')
grid on;

figure(2)
plot(time_line_landing, run_mean_queue_len_landing)
hold on;
yline(analytic_ss_queue_len, '--r');
hold off;
legend('MC Simulation Lq (Running Mean)', 'Analytic SS Lq');
xlabel('Time (minutes)')
ylabel('Mean Queue Length for Landing Planes')
title('Mean Queue Length Over Time for Landing Planes')
grid on;

figure(3)
plot(cum_arr_time_landing, q_wait_time_landing)
hold on;
yline(analytic_ss_queue_wait_time, '--r');
yline(mean_q_wait_time_landing, '--g');
hold off;
legend('MC Simulation Wq(i)', 'Analytic SS Wq', 'MC Simulation Wq (Mean)');
xlabel('Arrival time for Landing Planes')
ylabel('Wait Time in Queue for Landing Planes')
title('Wait Time in Queue Upon Arrival for Landing Planes')
grid on;

figure(4)
plot(time_line_landing, queue_len_landing)
hold on;
yline(analytic_ss_queue_len, '--r');
yline(mean_q_length_landing, '--g');
hold off;
legend('MC Simulation Lq(i)', 'Analytic SS Lq', 'MC Simulation Lq (Mean)');
xlabel('Time (minutes)')
ylabel('Queue Length for Landing Planes')
title('Queue Length Over Time for Landing Planes')
grid on;

figure(5)
plot(cum_arr_time_landing, q_wait_time_landing)
hold on;
yline(analytic_ss_queue_wait_time, '--r');
yline(mean_q_wait_time_landing, '--g');
hold off;
legend('MC Simulation Wq(i)', 'Analytic SS Wq', 'MC Simulation Wq (Mean)');
xlabel('Arrival time for Landing Planes')
ylabel('Wait Time in Queue for Landing Planes')
title('Wait Time in Queue Upon Arrival for Landing Planes')
grid on;

%--takeoff
figure(6)
plot(cum_arr_time_takeoff, run_mean_q_wait_time_takeoff)
hold on;
yline(analytic_ss_queue_wait_time, '--r');
hold off;
legend('MC Simulation Wq', 'Analytic SS Wq');
xlabel('Arrival time for Take-Off Planes')
ylabel('Mean Time in Queue for Take-Off Planes')
title('Mean Time in Queue Upon Arrival for Take-Off Planes')
grid on;

figure(7)
plot(time_line_takeoff, run_mean_queue_len_takeoff)
hold on;
yline(analytic_ss_queue_len, '--r');
hold off;
legend('MC Simulation Lq (Running Mean)', 'Analytic SS Lq');
xlabel('Time (minutes)')
ylabel('Mean Queue Length for Take-Off Planes')
title('Mean Queue Length Over Time for Take-Off Planes')
grid on;

figure(8)
plot(cum_arr_time_takeoff, q_wait_time_takeoff)
hold on;
yline(analytic_ss_queue_wait_time, '--r');
yline(mean_q_wait_time_takeoff, '--g');
hold off;
legend('MC Simulation Wq(i)', 'Analytic SS Wq', 'MC Simulation Wq (Mean)');
xlabel('Arrival time for Take-Off Planes')
ylabel('Wait Time in Queue for Take-Off Planes')
title('Wait Time in Queue Upon Arrival')
grid on;

figure(9)
plot(time_line_takeoff, queue_len_takeoff)
hold on;
yline(analytic_ss_queue_len, '--r');
yline(mean_q_length_takeoff, '--g');
hold off;
legend('MC Simulation Lq(i)', 'Analytic SS Lq', 'MC Simulation Lq (Mean)');
xlabel('Time (minutes)')
ylabel('Queue Length for Take-Off Planes')
title('Queue Length Over Time for Take-Off Planes')
grid on;

figure(10)
plot(cum_arr_time_takeoff, q_wait_time_takeoff)
hold on;
yline(analytic_ss_queue_wait_time, '--r');
yline(mean_q_wait_time_takeoff, '--g');
hold off;
legend('MC Simulation Wq(i)', 'Analytic SS Wq', 'MC Simulation Wq (Mean)');
xlabel('Arrival time for Take-Off Planes')
ylabel('Wait Time in Queue for Take-Off Planes')
title('Wait Time in Queue Upon Arrival for Take-Off Planes')
grid on;
