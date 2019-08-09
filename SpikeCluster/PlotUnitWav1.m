function PlotUnitWav(Handles)

Unit = Handles.Unit;
T = (0:size(Handles.Wav, 2)-1)/Handles.Para.Fs*1000;
XLim = [min(T)-T(2) max(T)+T(2)];
MaxWav = max(Handles.Wav, [], 2);
MinWav = min(Handles.Wav, [], 2);
YLim = [mean(MinWav)-7*std(single(MinWav))  mean(MaxWav)+7*std(single(MaxWav))];
% YLim = single([min(min(Handles.Wav))-10 max(max(Handles.Wav))+10]);
Unit(Unit<0) = 0;
Text = {'Unknown';'Unit1';'Unit2';'Unit3';};
for i = 0:3  % 0/1/2/3/4 -- Unknown/1/2/3...
    set(Handles.(['Unit' num2str(i)]), 'String', '');
    HCA = Handles.(['Wav' num2str(i) 'Axes']);
    cla(HCA, 'reset');
    set(gcf, 'CurrentAxes', HCA);
    switch get(Handles.WavPop, 'Value')
        case 1
            CurrWav = Handles.Wav(Unit==i,:);
        case {2;3;4;}
            CurrWav = Handles.Wav((Unit==i)&(~Handles.Outlier),:);
    end    
    if ~isempty(CurrWav)
        set(Handles.(['Unit' num2str(i)]), 'String', [Text{i+1} ' #' num2str(size(CurrWav, 1))]);
        if i==3
            NumHidden = sum(Unit>3);
            set(Handles.(Text{i+1}), 'String', [Text{i+1} ' #' num2str(size(CurrWav, 1)) ' #' num2str(NumHidden)]);
        end
        switch get(Handles.WavPop, 'Value')
            case 1  % Raw waveforms+Mean waveforms
                plot(T, CurrWav', Handles.Clr(i+1), 'LineWidth', .5);
                hold on;
                % Mean waveforms
                CurrMu = mean(CurrWav, 1);
                CurrSD = std(CurrWav, 0, 1);
                plot(T, CurrMu, Handles.Clr(2-min(1, i)), 'LineWidth', 2);
                hold on;
                plot(T, [CurrMu-CurrSD; CurrMu+CurrSD]', Handles.Clr(2-min(1, i)), 'LineWidth', .5);
                axis([XLim YLim]);
            case {2;3;4;}  % Mean waveforms
                CurrMu = mean(CurrWav, 1);
                CurrSD = std(CurrWav, 0, 1);
                hold on;
                %patch([T fliplr(T)], [max(CurrWav, [], 1) fliplr(min(CurrWav, [], 1))], Handles.ShadingClr{i+1}, 'EdgeColor', 'none');
                WavDensity(CurrWav, 1/30, Handles.Fs, YLim, Handles.Clr(i+1));
                %plot(T, CurrMu, Handles.Clr(2-min(1, i)), 'LineWidth', 2);
                %plot(T, [CurrMu-CurrSD; CurrMu+CurrSD]', Handles.Clr(2-min(1, i)), 'LineWidth', .5);
                axis([XLim YLim]);
                box off;
        end
    end
end