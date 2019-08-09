function PlotContSignal(Handles)

HCA = Handles.ContSignalAxes;
cla(HCA, 'reset');
set(gcf, 'CurrentAxes', HCA);
Fs = Handles.Fs;
StartT = Handles.StartT;
EndT = Handles.EndT;
% TTLT = Handles.TTLT;
TrialID = Handles.TrialID;
Thre = Handles.Thre(TrialID,:);
DetPolar = Handles.Para.DetPolar;

XLim = []; AbsT = [];
% XLim = []; AbsT = []; RelT = []; StimT = [];
for i = 1:numel(TrialID)
    switch get(Handles.ContSignalPop, 'Value');
        case 1  %Wideband Signal
            Signal = Handles.WBSignal{TrialID(i)};
        case 2  %Spike Signal
            Signal = Handles.SpikeSignal{TrialID(i)};
        case 3  %Transformed Spike Signal
            Signal = Handles.TSpikeSignal{TrialID(i)};
        case 4  %LFP
            if ~isfield(Handles, 'LFPSignal')
                return;
            end
            Signal = Handles.LFPSignal{TrialID(i)};
    end
    if isempty(Signal)
        continue;
    end
    T = (StartT(TrialID(i)):1/Fs:EndT(TrialID(i)));  %in sec
%     T = (StartT(TrialID(i)):1/Fs:EndT(TrialID(i))) +TTLT(TrialID(i),1);  %in sec
    ColID = 1:min(numel(T), numel(Signal));
    T = T(ColID); Signal = Signal(ColID);
    AbsT = [AbsT T];
%     RelT = [RelT T-TTLT(TrialID(i),1)];
%     StimT = [StimT TTLT(TrialID(i),1)];
    plot(T, Signal, 'Color', [.5 .5 .5], 'LineWidth', .5);
    hold on;
    CurrXLim = [XLim; [min(T) max(T)]];
    XLim(1,1) = min(CurrXLim(:,1));
    XLim(1,2) = max(CurrXLim(:,2));
end
xlim(XLim);
if get(Handles.ContSignalPop, 'Value')==3
    XLim = get(gca, 'xlim');
    for i = 1:numel(Thre)
        if DetPolar==3
            Clr = {[0 0 0];[.5 .5 .5]};
        else
            Clr = {[0 0 0];[0 0 0]};
        end
        plot(XLim, [Thre(i) Thre(i)], 'Color', Clr{i});
    end
end
% 
% % Mark the trial onset time
% YTick = get(gca, 'ytick')';
% plot([StimT; StimT], repmat(YTick([1 end]), [1 numel(StimT)]), 'Color', [.5 .5 .5], 'LineStyle', ':');
% 
% % Use relative time within trial as ticks
% % XTick = get(gca, 'xtick')';
% % [~, ID] = min(abs(repmat(AbsT, [numel(XTick) 1])-repmat(XTick, [1 numel(AbsT)])), [], 2);
% % set(gca, 'xtick', XTick, 'xticklabel', RelT(ID));