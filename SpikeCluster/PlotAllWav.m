function PlotAllWav(Handles)

switch get(Handles.WavPop, 'Value')
    case {1;2;}  % Waveforms
        Unit = Handles.Unit;
        Unit(Unit<0) = 0;
        UnitPool = unique(Unit);
    case {3;4;}  % Coefficients
        Unit = Handles.Unit(Handles.WavNum);
        % Unit(Unit<0) = 0;
        UnitPool = unique(Unit);
        UnitPool(UnitPool<0) = [];
end
T = (0:size(Handles.Wav, 2)-1)/Handles.Para.Fs*1000;
XLim = [min(T)-T(2) max(T)+T(2)];
MaxWav = single(max(Handles.Wav, [], 2));
MinWav = single(min(Handles.Wav, [], 2));
YLim = [median(MinWav)-7*std(MinWav)  median(MaxWav)+7*std(MaxWav)];
% YLim = single([min(min(Handles.Wav))-10 max(max(Handles.Wav))+10]);
cla(Handles.WavAxes, 'reset');
set(gcf, 'CurrentAxes', Handles.WavAxes);
switch get(Handles.WavPop, 'Value')
    case 1  % Raw waveforms
        for i = 1:numel(UnitPool)
            CurrWav = Handles.Wav(Unit==UnitPool(i),:);
            if ~isempty(CurrWav)
                ClrIdx = mod(UnitPool(i)-1, 5)+2;
                if UnitPool(i)==0
                    ClrIdx = 1;
                end
                plot(T, CurrWav', Handles.Clr(ClrIdx), 'LineWidth', .5);
                hold on;
                axis([XLim YLim]);
            end
        end
    case 2  % Mean waveforms
        for i = 1:numel(UnitPool)
            CurrWav = Handles.Wav((Unit==UnitPool(i))&(~Handles.Outlier),:);
            if ~isempty(CurrWav)
                ClrIdx = mod(UnitPool(i)-1, 5)+2;
                if UnitPool(i)==0
                    ClrIdx = 1;
                end
                CurrMu = mean(CurrWav, 1);
                CurrSD = std(CurrWav, 0, 1);
                hold on;
                patch([T fliplr(T)], [max(CurrWav, [], 1) fliplr(min(CurrWav, [], 1))], Handles.ShadingClr{ClrIdx}, 'EdgeColor', 'none');
                plot(T, CurrMu, Handles.Clr(ClrIdx), 'LineWidth', 2);
                plot(T, [CurrMu-CurrSD; CurrMu+CurrSD]', Handles.Clr(ClrIdx), 'LineWidth', .5);
                axis([XLim YLim]);
            end
        end
    case 3  % 2-D projection
        for i = 1:numel(UnitPool)
            CurrCoeff = Handles.Coeff(Unit==UnitPool(i),1:2);
            ClrIdx = mod(UnitPool(i)-1, 5)+2;
            if UnitPool(i)==0
                ClrIdx = 1;
            end
            scatter(CurrCoeff(:,1), CurrCoeff(:,2), 10, Handles.Clr(ClrIdx), 'fill', 'o', 'MarkerEdgeColor', 'w');
            hold on;
        end
    case 4  % 3-D projection
        if size(Handles.Coeff, 2)>=3
            for i = 1:numel(UnitPool)
                CurrCoeff = Handles.Coeff(Unit==UnitPool(i),1:3);
                ClrIdx = mod(UnitPool(i)-1, 5)+2;
                if UnitPool(i)==0
                    ClrIdx = 1;
                end
                scatter3(CurrCoeff(:,1), CurrCoeff(:,2),CurrCoeff(:,3), 10, Handles.Clr(ClrIdx), 'fill', 'o', 'MarkerEdgeColor', 'w');
                hold on;
            end
        end
end