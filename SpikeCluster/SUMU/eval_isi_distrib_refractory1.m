function [is_isi_su, perc_low_3ms] = eval_isi_distrib_refractory(...
    cluster_class, cur_spike_inds)
% eval_isi_distrib    
%

% Author: Ariel Tankus.
% Created: 20.12.2008.
% Modified: 25.05.2010.  V1.
% Copyright (C) Ariel Tankus, 2010.  All rights reserved.


perc_low_3ms_thresh = 3;   % percents of spikes allowed within 3ms after
                           % the previous spike (refractory period) [%].
                           % If the threshold is exceeded, the unit is a
                           % multiunit.

times = diff(cluster_class(cur_spike_inds, 2));

hist_edges = 0:100;
if (~isempty(times))    % empty times happens for a 1 spike cluster.
    N = histc(times, hist_edges);
else
    N = zeros(1, length(hist_edges));
end
% ARIEL: Modified, 02.11.2005:
aux_num = length(cur_spike_inds);
sum_low_3ms  = sum(N(1:3));
perc_low_3ms = sum_low_3ms ./ aux_num .* 100;

is_isi_su = (perc_low_3ms <= perc_low_3ms_thresh);
