function [susp_su, main_rise_inds, is_sim, avg_std_area_per_rise_len, ...
          rel_max_to_min] = is_susp_su(mean_waveform, std_waveform, ...
            var_waveform, par, with_plot)
% is_susp_su    Test whether a given waveform (assumed to be a mean
%               waveform of a class, NOT a single raw spike) should
%               be suspected as a single unit.
%
%               [susp_su, main_rise_inds, is_sim] = is_susp_su(mean_waveform,...
%                                                             std_waveform, par)
%               mean_waveform - 1xn - mean waveform of the spikes
%                                     in a single (assumed) class.
%               std_waveform  - 1xn - std. dev. of the waveform of the spikes
%                                     in a single (assumed) class.
%               par           - struct - classification parameters
%                                     stored in times_CSC files.
%                                     Should contain the w_pre field,
%                                     so we know at which point the
%                                     spikes are aligned (``peaked'').
%               susp_su - 1x1 - logical - true iff the waveform should be
%                                     suspected as a single unit
%                                     according to its waveforms only.
%                                     No attempt is made to evaluate the ISI.
%               main_rise_inds - 1x2 - [start, end] - indices of
%                                     the main rise in voltage,
%                                     starting the spike.
%               is_sim - 1x1 - logical - true iff the mean_waveform comes
%                                     from a simulation, not a real recording
%                                     (i.e., its values are in [-1, 1]).
%
%               See also: eval_class_quality, reg_times_csc,
%                         detect_main_voltage_rise.

% Author: Ariel Tankus.
% Created: 02.09.2008.
% Copyright (C) Ariel Tankus, 2010.  All rights reserved.


%rise_area_per_len_th = 0.09;   % if average area per rise length is below this
rise_area_per_len_th = 3;  % if average area per rise length is below this
                              % threshold, we suspect it is a single
                              % unit (we still have to check the ISI).
rel_max_to_min_th = 1.5;      % if the ratio of the height of the main rise
                              % relative to the beginning of the main rise
                              % and the minimum after the rise relative to the
                              % beginning of the main rise is below this
                              % threshold, the waveform shows small
                              % resemblence to spike, and is therefore
                              % declared multiunit.


% detect and mark the main rise:
if (mean_waveform(par.w_pre) > 0)
    if (~with_plot)
        [main_rise_inds, is_sim] = detect_main_voltage_rise(mean_waveform, par);
    else
        [main_rise_inds, is_sim, ind_bounds] = ...
            detect_main_voltage_rise(mean_waveform, par);
    end
else
    % invert waveform for negative spikes:
    if (~with_plot)
        [main_rise_inds, is_sim] = detect_main_voltage_rise(-mean_waveform, par);
    else
        [main_rise_inds, is_sim, ind_bounds] = ...
            detect_main_voltage_rise(-mean_waveform, par);
    end
end
if (with_plot)
    create_main_rise_detection_fig_seq;
    keyboard
end

if (isnan(main_rise_inds(1)))
    % suspected not-a-spike:
    susp_su = false;
    avg_std_area_per_rise_len = -1;    % smaller than any positive threshold.
    rel_max_to_min = 0;
    return;
end

if (is_sim)
    mean_waveform = mean_waveform .* 100;
end
main_rise = mean_waveform(main_rise_inds(1):main_rise_inds(2));

% compute the average of the area bounded between +-1std per unit length
% of the main rise:
std_area = 2.*sum(std_waveform(main_rise_inds(1):main_rise_inds(2)));
%%%var_area = sqrt(sum(var_waveform(main_rise_inds(1):main_rise_inds(2))));
%main_rise_vlen = sum(abs(diff(main_rise)));      % length of the rise.
main_rise_vlen = abs(main_rise(end)-main_rise(1));      % length of the rise.
if (main_rise_vlen == 0)
    fprintf('main_rise_vlen == 0\n');
%    keyboard
end
avg_std_area_per_rise_len = std_area ./ main_rise_vlen; 
%%%avg_std_area_per_rise_len = var_area ./ main_rise_vlen;
susp_su = (avg_std_area_per_rise_len);% < rise_area_per_len_th);
%if (~issorted(main_rise))
%    keyboard
%end

% Find minimal voltage after the rise:
after_rise = mean_waveform((main_rise_inds(2)+1):end);
min_after_rise = min(after_rise);

if (min_after_rise == main_rise(1))
    rel_max_to_min = log2(1000);           % very high.
else
    rel_max_to_min = log2(abs((main_rise(end) - main_rise(1)) ./ ...
                              (min_after_rise - main_rise(1))));
end

%if (susp_su)
%    if (rel_max_to_min < rel_max_to_min_th)
%        susp_su = false;
%    end
%end
