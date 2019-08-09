function PlotModel(Handles)
% Plot the model for clustering

if ~(size(Handles.Wav, 1)>Handles.MinNumWav)
    return;
end
cla(Handles.IndexAxes, 'reset');
set(gcf, 'CurrentAxes', Handles.IndexAxes);
switch Handles.Para.ClusterMethod
    case {1;5;}  % Single-linkage/SPC, plot tree
        [~, Col] = size(Handles.Model);
        for i = 1:min(Col-4, 12)
            ClrIdx = mod(i-1, 5)+2;
            semilogy(Handles.Model(:,2), Handles.Model(:,i+4), Handles.Clr(ClrIdx), 'LineWidth', 1);
            hold on;
        end
        XLim = [0 max(Handles.Model(:,2))];
        YLim = get(gca, 'YLim');
        semilogy(XLim, [Handles.MinSize Handles.MinSize], Handles.Clr(1), 'LineStyle', ':', 'LineWidth', 1);
        hold on;
        Temp = Handles.Model(Handles.Temp,2);
        semilogy([Temp Temp], YLim, Handles.Clr(1), 'LineStyle', ':', 'LineWidth', 1);
        hold on;
        axis([XLim YLim]);
    otherwise  % K-means/EM/tdistEM/AutomaticEM, plot 2-D projections
        Unit = Handles.Unit(Handles.WavNum);
%         Unit(Unit<0) = 0;
        UnitPool = unique(Unit);
        UnitPool(UnitPool<0) = [];
        for i = 1:numel(UnitPool)
            CurrCoeff = Handles.Coeff(Unit==UnitPool(i),1:2);
            ClrIdx = mod(UnitPool(i)-1, 5)+2;
            if UnitPool(i)==0
                ClrIdx = 1;
            end
            scatter(CurrCoeff(:,1), CurrCoeff(:,2), 10, Handles.Clr(ClrIdx), 'fill', 'o');
            hold on;
        end
end