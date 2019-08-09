Copyright (C) Ariel Tankus, 2010.  All rights reserved.

The algorithm of SUMU has been described in detail in:

@Article{ariel:sumu,
  author =       {Ariel Tankus and Yehezkel Yeshurun and Itzhak Fried},
  title =        {An automatic measure for classifying clusters of suspected
                  spikes into single cells versus multiunits},
  journal =      {Journal of Neural Engineering},
  year =         {2009},
  volume =       {6},
  number =       {5},
  month =        {Oct.},
  pages =        {056001},
}

If you use the SUMU algorithm and publish the results of its run (i.e., the
classification of clusters into single units vs. multiunits) in any form,
you are REQUIRED to cite the aforementioned paper in your publication.
This includes the case where you publish any statistic of the classification
or any other derivative form of the classification (for example, if you say:
"We had 400 single units and 500 multiunits" and their classification was
obtained by the SUMU algorithm, then you MUST cite the aforementioned paper).
The citation MUST appear in EVERY publication which uses the classification
results of SUMU, whether it is a direct or an indirect usage (i.e., citation of
the aforementioned paper in one of your publications and then citing your own
paper and saying something like: "We used the same classification method as in
our previous paper" is UNACCEPTABLE.  In this example, the second paper violates
the copyright agreement.  The aforementioned paper should be cited in EVERY
publication.  In this example both the first and second papers should cite the
aforementioned paper.)

The code is provided AS IS, without any warranty.  Your using the code in this
distribution is done AT YOUR OWN RISK.  Copying, modification and distribution
of the code are strictly prohibited, except for personal usage.  Any
modification of the code for personal use is subject to the same copyright terms
of the original distribution (i.e., if you publish something that was created by
running the modified code, you still need to cite the aforementioned paper with
exactly the same terms as if you ran the original code itself only).

If you cannot fulfil the terms of this license for any reason, you are required
to delete any copy of the software and not use it in any way.


Running the code:
-----------------

The SUMU algorithm was created for files in the format of WaveClus (a
spike sorting method by Quian-Quiroga et al.).  For files that were sorted by
WaveClus, running the SUMU algorithm is as simple as running:

eval_class_quality_session(ch_list);

where ch_list is a 1xn list of channel IDs (e.g., ch_list = 1:64; ).

---
For data not sorted with WaveClus, run:

[is_clear_spike, is_spike, spike_inds_array, params_values] = ...
        eval_class_quality(cluster_class, par, spikes, with_plot);

where:
cluster_class - nx2 - [cluster_ID, spike_time].  The cluster ID is within the
                      current channel, where 0 means not-a-spike (i.e., trash)
                      and 1:k are the clusters decided upon for the current
                      channel.
par - struct - should contain the field: par.w_pre .  This field is a 1x1 scalar
                      denoting the number of sampling points at which the spike
                      is aligned (i.e., the peak of the spike).
spikes - nxl - is the spike waveform of each spike in cluster_class.  Each
                      waveform contains l samples (note that par.w_pre < l).
with_plot - 1x1 - logical - true if you want plots in the process; false,
                      otherwise.
is_clear_spike - rx1 - logical - true: The cluster is a single unit.
                      false: The cluster is a multiunit.  (r: #clusters,
                      excluding the trash.)
is_spike - rx1 - logical - true: The cluster represents a cluster of spikes.
                      false: The cluster does NOT represent a cluster of spikes,
                      and should be added to the "trash" cluster.  A suspected
                      cluster is denoted as "false" only in extreme cases where
                      its shape is extremely different from a typical spike and
                      the main rise in voltage can not be detected.
spike_inds_array - rx1 - cell - each cell contains the indices of the spikes in
                      the corresponding cluster (indices into cluster_class).
params_values - rx2 - [avg_std_area_per_rise_len, perc_low_3ms] -
                      The first column is the average area of the standard
                      deviation around the main rise per rise length (see the
                      papers for details about this feature).
                      The second column is the percentage of inter-spike
                      intervals that are shorter than 3ms.

Changing the threshold:
-----------------------
In is_susp_su.m the threshold on the relative area around the main rise in
voltage is set (the first code line).  The default is:
rise_area_per_len_th = 3;
but you may change it for your own needs.

---

-- Ariel.
25.05.2010.
