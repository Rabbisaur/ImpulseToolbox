function [Output] = ExtractFeature(Input, Method, varargin)
%% Check INPUT
if ~isnumeric(Input) || size(Input, 2)<2
    error('INPUT must be a numeric matrix with 2/more coloumns.');
end
WaveletType = lower(varargin{1});
WaveletScale = varargin{2};
%% Extract discriminative features
NumRow = size(Input, 1);
switch Method
    case 1  %raw waveform
        Output = Input;
    case 2  %reduced feature set, global maximum/minimum/width
        Output = nan([NumRow 3], 'single');
        [Output(:,1) MaxIndex] = max(Input, [], 2);
        [Output(:,2) MinIndex] = min(Input, [], 2);
        Output(:,3) = MaxIndex-MinIndex;
    case 3  %PCA
        Input = single(Input);
        [~, Output] = princomp(Input);
    case 4  %wavelet transform(Quian Quiroga, Neural Computation2004)
        for i = NumRow:-1:1
            Output(i,:) = wavedec(single(Input(i,:)), WaveletScale, WaveletType);
        end
    case 5  %wavelet packets decomposition(Hulata, J Neurosci. Methods2002)
        Input = single(Input);
        Output = Input;
        WPTree = wpdec(mean(Input,1), WaveletScale, WaveletType);
        WPTree = besttree(WPTree);
        WPNode = get(WPTree, 'tn');
        for i= 1:NumRow
            WPTree = wpdec(Input(i,:), WaveletScale, WaveletType);
            NumCoeff = 0;
            for j = 1:numel(WPNode)
                Coeff = wpcoef(WPTree, WPNode(j));
                Output(i,NumCoeff+1:NumCoeff+numel(Coeff)) = Coeff;
                NumCoeff = NumCoeff+numel(Coeff);
            end
        end
    case 6  %multiwavelets transform(Geng, Neurocomputing2010)
        ScalingCoeff = floor(size(Input, 2)/2^WaveletScale);
        for i=NumRow:-1:1
            PreInput = prep1D_rr(Input(i,:), WaveletType);
%             PreInput = prep1D_appe(Input(i,:), [WaveletType 'ap']);
            PreInput = single(PreInput);
            Coeff = dec1D_pe(PreInput, WaveletType, WaveletScale);
            Coeff = Coeff(:,ScalingCoeff+1:end);  %multiwavelets coeff
            Output(i,:)=Coeff(:);
        end
    otherwise
        error('No current Method was found.');
end