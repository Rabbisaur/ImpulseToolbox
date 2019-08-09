function LLH = LLHWeibullComputeFit(FitP,Lum,nTrials,Seen,nCatch,nCR)
  pFA = FitP(1);
  MaxPerf = FitP(2);
  Slope = FitP(3);
  Thresh = FitP(4);
  Perf = MaxPerf + (pFA-MaxPerf) * exp(-(Lum/Thresh).^Slope);
  NLums = size(Lum);
  NLums = NLums(2);
  nMiss = nTrials - Seen;
  LLH2=0;
  nFA=nCatch-nCR;
  for Lum = 1:NLums
      PCorr = Perf(Lum);
      NTr = nTrials(Lum);
      Seen1 = Seen(Lum);
      Miss1 = nMiss(Lum);
      LLH1 = log(PCorr)*Seen1 + log(1-PCorr)*Miss1; % + log(nchoosek(NTr,Seen1))
    LLH2 = LLH2 - LLH1;     
  end
  PCorr = pFA;
  LLH1 = log(PCorr)*nFA + log(1-PCorr)*nCR; % + log(nchoosek(nCatch,nFA))
  LLH2 = LLH2 - LLH1;  
  LLH = LLH2;
end
