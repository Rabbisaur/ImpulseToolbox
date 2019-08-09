function Para = DefaultPara() 
% Set default parameters
Para.Load = 1;
Para.ElecNum = 1;
Para.Fs = 30000;
Para.SegDuration = 60;
Para.Filter = 0;
Para.HipassFreq = 300;
Para.LopassFreq = 8000;
Para.FilterDirection = 1;
Para.HipassType = 1;
Para.HipassOder = 4;
Para.LopassType = 1;
Para.LopassOrder = 3;
Para.Det = 0;
Para.DetMethod = 3;
Para.SDMethod = 5;
Para.SDDuration = 5;
Para.DetThre = 2;
Para.DetPolar = 3;
Para.SpikeDuration = 1.6;
Para.RefractoryPeriod = 3;
Para.Align = 0;
Para.ReverseFilterWavPop = 1;
Para.InterpMethod = 1;
Para.InterpReso = 15;
Para.AlignMethod = 2;
Para.AlignTime = 1/3;
Para.Extract = 1;
Para.NumClusteredSpike = 5000;
Para.ExtractMethod = 6;
Para.WaveletType = 1;
Para.WaveletScale = 4;
Para.Select = 1;
Para.SelectMethod = 3;
Para.VarRatio = 90;
Para.NumCoeff = 10;
Para.Cluster = 1;
Para.OutlierRemoval = 2;
Para.ClusterMethod = 5;
Para.NumCluster = 2;
Para.MinClusterRatio = 2.5;
Para.AutoTemplate = 2;
Para.AutotempSNRThre = 2;
Para.AssignRemainder = 2;
Para.AssignOutlier = 0;
Para.HMM = 0;
Para.SNR = 1;
Para.IsolationScore = 0;
Para.FalseNeg = 0;
Para.FalsePos = 0;
Para.SingleMultiUnit = 1;
Para.FalseNegDet = 1;
Para.FalsePosRP = 1;
Para.FalseNegCluster = 1;
Para.FalsePosCluster = 1;
Para.FalseNegCensor = 1;
Para.LRatio = 1;
Para.IsolationDist = 1;
Para.IsolationInfo = 1;
% Pop options
Para.FilterDirectionOpt = {'Forward & Reverse';'Forward Only';'Reverse Only'};
Para.HipassOpt = {'None';'Butterworth';'Chebyshev Type I';'Chebyshev Type II';'Elliptic'};
Para.LopassOpt = {'None';'Butterworth';'Chebyshev Type I';'Chebyshev Type II';'Elliptic'};
Para.DetMethodOpt = {'Raw Voltage';'Nonlinear Energy Operator';'Mathematical Morphonogy';};
Para.SDMethodOpt = {'Root mean square';'Median absolute deviation';'Cap fitting';'Duty-cycle keeping';'Max-min spread';};
Para.DetPolarOpt = {'Negative';'Positive';'Either';'Both';};
Para.ReverseFilterWavOpt = {'No';'Yes'};
Para.InterpMethodOpt = {'Raw Waveform';'Cubic Spline';'MATLAB V5 Cubic';'Linear Spline';'Nearest Neighbor';};
Para.AlignMethodOpt = {'Raw Waveform';'Global Minimum';'Global Maximum';'Global Absolute Maximum';'Global Minimum Slope';...
                       'Global Maximum Slope';'Global Absolute Slope';'Multi-peak Energy Comparison';};
Para.ExtractMethodOpt = {'Raw Waveform';'Reduced Feature Set';'Principal Component Analysis';'Wavelet Coefficients';...
                         'Wavelet Packets Coefficients';'Multiwavelets Coefficients';};
Para.WaveletTypeOpt = {'Haar';'Db2';'Bior1.1';'Coif5';
                       'Db3';'Db4';'Db5';'Db10';'Db20';'Db30';'Db40';'Coif1';'Coif1';'Coif2';'Coif3';'Coif4';...
                       'Sym4';'Sym5';'Sym10';'Sym20';'Sym30';'Sym40';'Dmey';'Bior1.3';'Bior1.5';'Bior2.2';...
                       'Bior2.4';'Bior2.6';'Bior2.8';'Bior3.1';'Bbior3.3';'Bbior3.5';'Bbior3.7';'Bior3.9';...
                       'Bbior4.4';'Bbior5.5';'Bbior6.8';'Rbio1.1';'Rbio1.3';'Rbio1.5';'Rbio2.8';'Rbio3.1';...
                       'Rbio3.3';'Rbio3.5';'Rbio3.7';'Rbio3.9';'Rbio4.4';'Rbio5.5';'Rbio6.8';};
Para.MultiwaveletTypeOpt = {'Bighm2';'Bih34n';'Bighm6';'La8';'Cardbal4';...
                            'Haar';'D4';'Bi9';'Bi7';'Bi5';'Bi3';'GHM';'CL';'Sa4';'Bih52s';'Bih32s';'Bih54n';...
                            'Cardbal2';'Cardbal3';'Clbal';'GHMbal';};
Para.SelectMethodOpt = {'All Coefficients';'Maximum Variance';'KS Statistics for Normality';'Square Sum of eCDF-GaussianCDF difference';'GOF for Gaussian Mixtures';};
Para.OutlierRemovalOpt = {'None';'Global Outliers';'Outliers in Separate Clusters';'Global & Cluster Outliers';};
Para.ClusterMethodOpt = {'Single Linkage';'K-Means';'Expectation-maximization';'t-Distribution Expectation-maximization';'Superparamagnetic';'Automatic EM';};
Para.NumClusterOpt = {'Auto';'1';'2';'3';'4';'5';'6';'7';'8';};
Para.AutoTemplateOpt = {'Default Templates';'Peak-Valley Template';'Mathematical Morphonogy Template';'Peak-Valley Difference';'Mathematical Morphonogy';};
Para.AssignRemainderOpt = {'Never';'Waveform, nearest neigbor';'Waveform, center';'Waveform, maximum likehood';'Waveform, Mahalanobis';...
                           'Coefficient, nearest neigbor';'Coefficient, center';'Coefficient, maximum likehood';'Coefficient, Mahalanobis';};
Para.AssignOutlierOpt = {'Never';'Segment Matching';};
save('DefaultPara.mat', 'Para', '-v7.3');