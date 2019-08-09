function [is_clear_spike, params_values] = eval_class_quality_from_file(ch, ...
            times_prefix, with_plot)
% eval_class_quality_from_file    

% Author: Ariel Tankus.
% Created: 02.09.2008.
% Modified: 25.05.2010. V1.
% Copyright (C) Ariel Tankus, 2010.  All rights reserved.


if (nargin < 2)
    times_prefix = 'times';
end
if (nargin < 3)
    with_plot = false;
end

if (isnumeric(ch))
    filename = sprintf('%s_CSC%d.mat', times_prefix, ch);
else
    filename = sprintf('%s_%s.mat', times_prefix, ch);
end

if (~exist(filename, 'file'))
    is_clear_spike = [];
    params_values  = [];
    return;
end
load(filename);

[is_clear_spike, is_spike, spike_inds_array, params_values] = ...
    eval_class_quality(cluster_class, par, spikes, with_plot);
