function PlotProjections(Handles)

figure;
Pos = get(gcf, 'OuterPosition');
ScrSize = get(0, 'screensize');
Pos(1) = ScrSize(3)/2-Pos(3)/2;
Pos(2) = ScrSize(4)/2-Pos(4)/2;
set(gcf, 'Name', 'Coefficients Projections', 'OuterPosition', Pos);
set(gcf, 'outerposition', ScrSize);

Unit = Handles.Unit(Handles.WavNum);
% Unit(Unit<0) = 0;
UnitPool = unique(Unit);
UnitPool(UnitPool<0) = [];
NumCoeff = size(Handles.Coeff, 2);
for i = 1:NumCoeff
    for j = i+1:NumCoeff
        subplot(NumCoeff, NumCoeff, (i-1)*NumCoeff+j);
        hold on;
        for k = 1:numel(UnitPool)
            ClrIdx = mod(UnitPool(k)-1, 5)+2;
            if UnitPool(k)==0
                ClrIdx = 1;
            end
            scatter(Handles.Coeff(Unit==UnitPool(k),i), Handles.Coeff(Unit==UnitPool(k),j), ...
                20, Handles.Clr(ClrIdx), '.');
            axis off;
            axis square;
        end
    end
end