function Perf = LumWeibull(Lums, FA, MaxPerf, Slope, Thresh)
  Perf = MaxPerf + (FA-MaxPerf) * exp(-(Lums/Thresh).^Slope);
end

