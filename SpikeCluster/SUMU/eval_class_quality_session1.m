function eval_class_quality_session(ch_list)
% eval_class_quality_session    

% Author: Ariel Tankus.
% Created: 03.01.2009.
% Modified: 25.05.2010. V1.
% Copyright (C) Ariel Tankus, 2010.  All rights reserved.


if (nargin < 1)
    ch_list = 1:64;
end 

num_ch = length(ch_list);
sumu_per_ch = cell(num_ch, 1);
params_values_per_ch = cell(num_ch, 1);

for i=1:num_ch
    ch = ch_list(i);
    fprintf('Channel: %d\n', ch);

    [is_clear_spike, params_values] = eval_class_quality_from_file(ch, ...
                'times', false);
    if (isempty(is_clear_spike))
        % No clusters defined in a times_CSC file:
        continue;
    end
    
    sumu_per_ch{i} = is_clear_spike;
    params_values_per_ch{i} = params_values;
end

sumu_per_cl = cat(1, sumu_per_ch{:});
params_values_per_cl = cat(1, params_values_per_ch{:});

save sumu_is_su.mat sumu_per_ch sumu_per_cl params_values_per_ch ...
    params_values_per_cl;
