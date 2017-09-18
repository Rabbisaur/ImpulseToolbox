function PlotIndex(Handles)
% Plot statistics index for each unit

Unit = Handles.Unit;
Unit(Unit<0) = 0;
YLim = [0 0];
for i = 0:3  % 0/1/2/3/4 -- Unknown/1/2/3...
    HCA = Handles.(['Index' num2str(i) 'Axes']);
    cla(HCA, 'reset');
    set(gcf, 'CurrentAxes', HCA);
    set(Handles.(['Index' num2str(i)]), 'String', '');
    UnitIdx = Unit==i;
    switch get(Handles.IndexPop, 'Value')
        case 1  %Inter-spike Interval
            if sum(UnitIdx)>0
                ISI = zeros(size(Handles.Spike));
                ISI(UnitIdx) = [diff(Handles.Spike(UnitIdx)); 0]*1000;
                Idx = cumsum(Handles.NumSpikeTrial(Handles.NumSpikeTrial~=0));
                UnitIdx(Idx) = 0;
                CurrISI = ISI(UnitIdx);
                CurrISI = CurrISI(CurrISI>=0);
                if isempty(CurrISI)
                    continue;
                end
%                 hist(CurrISI(CurrISI<510), -10.5:1:510.5);
%                 HPatch = findobj(gca, 'Type', 'patch');
%                 set(HPatch, 'FaceColor', Handles.Clr(i+1), 'EdgeColor', Handles.Clr(i+1));
                Edge = logspace(log10(10^-10), log10(2000), 1000);                
                NumWav = histc(CurrISI, Edge);
                semilogx(Edge(1:end-1), NumWav(1:end-1), Handles.Clr(i+1), 'LineWidth', 1);
                hold on;
                YBound = get(gca, 'YLim');
                plot([3 3], YBound, [':' Handles.Clr(2-min(1, i))], 'LineWidth', 1);
                NumInRP = sum((CurrISI<2.5).*(CurrISI>=0));
                InRPRatio = NumInRP/numel(CurrISI)*100;
                set(Handles.(['Index' num2str(i)]), 'String', ['<2.5ms #'...
                    num2str(NumInRP) '--' num2str(InRPRatio, '%.2f') '%']);
                xlim([0 10000]);
            end
        case 2  %Detection/Threshold Histogram
            if sum(UnitIdx)>0
                CurrWav = single(Handles.Wav(UnitIdx,:));
                MinWav = min(CurrWav, [], 2);
                MuMin = mean(MinWav);
                SDMin = std(MinWav);
                Edge = linspace(MuMin-3*SDMin, MuMin+3*SDMin, 100);
                NumWav = histc(MinWav, Edge);
                plot(Edge(1:end-1), NumWav(1:end-1), Handles.Clr(i+1), 'LineWidth', 1);
                hold on;
                MaxWav = max(CurrWav, [], 2);
                MuMax = mean(MaxWav);
                SDMax = std(MaxWav);
                Edge = linspace(MuMax-3*SDMax, MuMax+3*SDMax, 100);
                NumWav = histc(MaxWav, Edge);
                plot(Edge(1:end-1), NumWav(1:end-1), [':' Handles.Clr(i+1)], 'LineWidth', 1);
                hold on;
            end
        case 3  %Amplitude Stability
            XLim = [0 max(max(Handles.Spike))];
            if sum(UnitIdx)>0
                CurrWav = single(Handles.Wav(UnitIdx,:));
                CurrSpike = Handles.Spike(UnitIdx,:);
                scatter(CurrSpike, min(CurrWav, [], 2), 20, Handles.Clr(i+1), '.');
                hold on;
                scatter(CurrSpike, max(CurrWav, [], 2), 20, Handles.Clr(i+1), '+');
                hold on;
                YBound = get(gca, 'YLim');
                YLim(1) = min(YLim(1), YBound(1));
                YLim(2) = max(YLim(2), YBound(2));
            end
            xlim(XLim);
    end
end
for i = 0:3  % 0/1/2/3/4 -- Unknown/1/2/3...
    HCA = Handles.(['Index' num2str(i) 'Axes']);
    set(gcf, 'CurrentAxes', HCA);
    switch get(Handles.IndexPop, 'Value')
        case 3  %Amplitude Stability
            ylim(YLim);
    end
end