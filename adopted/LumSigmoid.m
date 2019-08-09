function Perf = LumSigmoid(Lums, FA, MaxPerf, Slope, Thresh)
  ELums = exp((Lums-Thresh)*Slope);
  Perf = FA + (MaxPerf-FA) * ELums./(1 + ELums);
end

