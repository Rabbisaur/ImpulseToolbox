%function [] = create_main_rise_detection_fig()
% create_main_rise_detection_fig    

% Author: Ariel Tankus.
% Created: 17.02.2009.
% Copyright (C) Ariel Tankus, 2010.  All rights reserved.


ind = ind_bounds(1);
ind_v_rise = ind_bounds(2);
max_ind = main_rise_inds(1);
inds_rise = main_rise_inds(1):main_rise_inds(2);

%t = -(21*33):33:(48*33);     % for samp. rate of 30303Hz
t = -(19*36):36:(44*36);     % for samp. rate of 27777Hz


m = mean_waveform;
m_p_s = mean_waveform + std_waveform;
m_m_s = mean_waveform - std_waveform;

t_range = [t(1), t(end)];
v_range = [min(m_m_s), max(m_p_s)];

t = 1:numel(m);
% Mean waveform:
figure;
h_mean = plot(t, m, '+-');
xlabel('Time [\mus]');
ylabel('Voltage [\muV]');
set(gca, 'XLim', t_range);
set(gca, 'YLim', v_range);

% Upper-bound:
figure;
t_inds = find(t <= 0);
h_mean = plot(t(t_inds), m(t_inds), '+-');
hold on;
h_upper = plot(t(ind_v_rise), m(ind_v_rise), 'r*');   % upper bound.
xlabel('Time [\mus]');
ylabel('Voltage [\muV]');
set(gca, 'XLim', t_range);
set(gca, 'YLim', v_range);

% Lower-bound:
figure;
h_mean = plot(t(1:ind_v_rise), m(1:ind_v_rise), '+-');
hold on;
h_upper = plot(t(ind_v_rise), m(ind_v_rise), 'r*');   % upper bound.
h_init = plot(t(ind), m(ind), 'g*');  % initial point est. 
xlabel('Time [\mus]');
ylabel('Voltage [\muV]');
set(gca, 'XLim', t_range);
set(gca, 'YLim', v_range);

% Curvature-based detection: 
figure;
h_mean = plot(t(ind:ind_v_rise), m(ind:ind_v_rise), '+-');
hold on;
h_upper = plot(t(ind_v_rise), m(ind_v_rise), 'r*');   % upper bound.
h_init  = plot(t(ind), m(ind), 'g*');  % initial point est. 
h_curv  = plot(t(max_ind), m(max_ind), 'k*');         % curvature-based.
xlabel('Time [\mus]');
ylabel('Voltage [\muV]');
set(gca, 'XLim', t_range);
set(gca, 'YLim', v_range);

% Complete draw of mean + std. dev. on main rise only:
figure;
h_sd = plot(t, m_p_s, '+-', 'Color', [0.7, 0.7, 0.7]);
hold on;
plot(t, m_m_s, '+-', 'Color', [0.7, 0.7, 0.7]);
line([t(inds_rise); t(inds_rise)], [m_m_s(inds_rise); m_p_s(inds_rise)], ...
     'Color', [0.7, 0.7, 0.7]);
h_mean = plot(t, m, '+-');
hold on;
h_upper = plot(t(ind_v_rise), m(ind_v_rise), 'r*');   % upper bound.
%fill_circle(ind_v_rise, mean_waveform(ind_v_rise), 0.7, [1, 0, 0]);% upper bound.
h_init = plot(t(ind), m(ind), 'g*');  % initial point est. 
%fill_circle(ind, mean_waveform(ind), 0.7, [0, 1, 0]);  % initial point est. 
h_curv = plot(t(max_ind), m(max_ind), 'k*');         % curvature-based.
%fill_circle(max_ind, mean_waveform(max_ind), 0.7, [0, 0, 0]);  % curvature-based.
legend([h_mean, h_sd, h_upper, h_init, h_curv], ...
       {'Mean Waveform', 'Mean \pm SD', 'Upper bound.', 'Lower bound.', 'Curvature-based.'});
xlabel('Time [\mus]');
ylabel('Voltage [\muV]');
set(gca, 'XLim', t_range);
set(gca, 'YLim', v_range);

% Enlargement of the main rise with std. dev. marked:
figure;
h_sd = plot(t(inds_rise), m_p_s(inds_rise), '-', 'Color', [0, 0, 0]);
hold on;
plot(t(inds_rise), m_m_s(inds_rise), '-', 'Color', [0, 0, 0]);
h_mean = plot(t(inds_rise), m(inds_rise), '-');
line([t(inds_rise); t(inds_rise)], [m_m_s(inds_rise); m_p_s(inds_rise)], ...
     'Color', [1, 0, 0]);
%h_curv = plot(t(max_ind), m(max_ind), 'k*');         % curvature-based.
xlabel('Time [\mus]');
ylabel('Voltage [\muV]');
set(gca, 'XLim', t_range);
set(gca, 'YLim', v_range);
%set(gca, 'XLim', [t(inds_rise(1)), t(inds_rise(end))]);
%set(gca, 'YLim', [min(m_m_s(inds_rise(1):inds_rise(end))), ...
%                  max(m_p_s(inds_rise(1):inds_rise(end)))]);

% Complete draw of mean + std. dev. on main rise only:
figure;
h_sd = plot(t, m_p_s, '-', 'Color', [0, 0, 0]);
hold on;
plot(t, m_m_s, '-', 'Color', [0, 0, 0]);
line([t(inds_rise); t(inds_rise)], [m_m_s(inds_rise); m_p_s(inds_rise)], ...
     'Color', [1, 0, 0]);
h_mean = plot(t, m, '-');
xlabel('Time [\mus]');
ylabel('Voltage [\muV]');
set(gca, 'XLim', t_range);
set(gca, 'YLim', v_range);
