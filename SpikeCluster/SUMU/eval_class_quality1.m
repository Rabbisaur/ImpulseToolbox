function [is_clear_spike, is_spike, spike_inds_array, params_values] = ...
        eval_class_quality(cluster_class, par, spikes, with_plot)
% eval_class_quality    Evaluate the quality of classification of a suspected
%                       cluster of spikes into single unit vs. multiunit
%                       according to the SUMU algorithm.
%
%                       [is_clear_spike, is_spike, spike_inds_array,
%                       params_values] = eval_class_quality(cluster_class, par,
%                                                           spikes, with_plot);
%                       cluster_class - nx2 - [cluster_ID, spike_time].  The
%                                             cluster ID is within the
%                                             current channel, where 0 means
%                                             not-a-spike (i.e., trash) and
%                                             1:k are the clusters decided
%                                             upon for the current channel.
%                       par - struct - should contain the field: par.w_pre .
%                                             This field is a 1x1 scalar
%                                             denoting the number of sampling
%                                             points at which the spike is
%                                             aligned (i.e., the peak of the
%                                             spike).
%                       spikes - nxl - is the spike waveform of each spike in
%                                             cluster_class.  Each waveform
%                                             contains l samples (note that
%                                             par.w_pre < l).
%                       with_plot - 1x1 - logical - true if you want plots in
%                                             the process; false, otherwise.
%                       is_clear_spike - rx1 - logical - true: The cluster is
%                                             a single unit.
%                                             false: The cluster is a
%                                             multiunit.  (r: #clusters,
%                                             excluding the trash.)
%                       is_spike - rx1 - logical - true: The cluster
%                                             represents a cluster of spikes. 
%                                             false: The cluster does NOT
%                                             represent a cluster of spikes,
%                                             and should be added to the
%                                             "trash" cluster.  A suspected
%                                             cluster is denoted as "false"
%                                             only in extreme cases where its
%                                             shape is extremely different
%                                             from a typical spike and the
%                                             main rise in voltage can not be
%                                             detected.
%                       spike_inds_array - rx1 - cell - each cell contains
%                                             the indices of the spikes in
%                                             the corresponding cluster
%                                             (indices into cluster_class).
%                       params_values - rx2 - [avg_std_area_per_rise_len,
%                                              perc_low_3ms] - 
%                                             The first column is the average
%                                             area of the standard deviation
%                                             around the main rise per rise
%                                             length (see the papers for
%                                             details about this feature).
%                                             The second column is the
%                                             percentage of inter-spike
%                                             intervals that are shorter than
%                                             3ms.
%                       See: 
% @Article{ariel:sumu,
%   author =       {Ariel Tankus and Yehezkel Yeshurun and Itzhak Fried},
%   title =        {An automatic measure for classifying clusters of suspected
%                   spikes into single cells versus multiunits},
%   journal =      {Journal of Neural Engineering},
%   year =         {2009},
%   volume =       {6},
%   number =       {5},
%   month =        {Oct.},
%   pages =        {056001},
% }
%
% for details of the algorithm.

% Author: Ariel Tankus.
% Created: 19.08.2008.
% Modified: 25.05.2010. V1.
% Copyright (C) Ariel Tankus, 2010.  All rights reserved.


if (nargin < 4)
    with_plot = false;
end

num_params = 2;

class_ids = unique(cluster_class(:, 1));
class_ids = setdiff(class_ids, 0);           % ignore trash class 0.
num_classes = length(class_ids);
waveform_len = size(spikes, 2);
spike_inds_array = cell(num_classes, 1);
mean_waveforms  = zeros(num_classes, waveform_len);
var_waveforms   = zeros(num_classes, waveform_len);
std_waveforms   = zeros(num_classes, waveform_len);
se_waveforms    = zeros(num_classes, waveform_len);
is_clear_spike  = zeros(num_classes, 1);
is_spike        = true(num_classes, 1);
params_values   = NaN(num_classes, num_params);

colors = {'r', 'g', 'b', 'y', 'c', 'k', 'm', 'r', 'g', 'b'};
legend_str = cell(num_classes, 1);
is_first_plot = true;

for i=1:num_classes

    % Basic cluster statistics:
    spike_inds_array{i} = find(cluster_class(:, 1) == class_ids(i));
    mean_waveforms(i, :) = mean(spikes(spike_inds_array{i}, :), 1);
    var_waveforms(i, :)  = var(spikes(spike_inds_array{i}, :), 0, 1);
    std_waveforms(i, :)  = std(spikes(spike_inds_array{i}, :), 0, 1);
    se_waveforms(i, :)   = std_waveforms(i, :) ./ ...
                               sqrt(length(spike_inds_array{i}));

    [susp_su, main_rise_inds, is_sim, avg_std_area_per_rise_len, ...
     rel_max_to_min] = is_susp_su(mean_waveforms(i, :), std_waveforms(i, :), ...
                                  var_waveforms(i, :), par, with_plot);
    is_clear_spike(i) = susp_su;
    is_spike(i) = (~isnan(main_rise_inds(1)));   % is waveform suspected
                                                 % as a spike, or can be
                                                 % rejected at this stage?

    % Initialize params to NaN, in case they are not evaluated due to
    % earlier params:
    perc_low_3ms = NaN;

    [is_isi_su, perc_low_3ms] = eval_isi_distrib_refractory(cluster_class, ...
                                                      spike_inds_array{i});
    %is_clear_spike(i) = (susp_su & is_isi_su);
    
    params_values(i, :) = [avg_std_area_per_rise_len, perc_low_3ms];

end
