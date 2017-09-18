function [elecRF,PSTH,FigurePath] = RFtuningElecFun(pathx,pathy,scancenter,barwidth,figurepath,elecunits,arraymap)

if ispc
    slash = '\';
else
    slash = '/';
end

sessionpath.x = pathx;
sessionpath.y = pathy;
unitfilename = 'unit.mat';

idx = find(sessionpath.x==slash,1,'last');
idx2 = find(sessionpath.x=='_',1,'last');
sessionname = [sessionpath.x(idx+1:idx2-1),sessionpath.x(idx2+1:end)];

% config
edges = ([0:10:1000] - 200)/1000; % unit is seconds.
time = ([0:10:1000] - 200)/1000;
time(end) = [];
timewindow = [0 0.5]; % stimulus on time = 0.5s
center = scancenter;
% barwidth = 0.3;
rsquarethres = 0.3;
xSigma = 1.96;


PSTH(1) = TSLgetPSTHelec(sessionpath.x,edges,elecunits,unitfilename);
% PSTH = TSLgetPSTHelec(sessionpath,edges,UnitSelect,unitfilename)
PSTH(2) = TSLgetPSTHelec(sessionpath.y,edges,elecunits,unitfilename);
% PSTH(3) = TSLgetPSTHelec(sessionpath.o,edges,elecunits);

g = fittype(@(a,b,c,d,x) a.*exp(-((x-b).^2/(2.*c.^2)))+d,'independent','x',...
    'dependent','y');

for i = 1:2 % x, y, o
    NumConditions = PSTH(i).NumConditions;
    NumElec = PSTH(i).NumElec;
    if i < 3
        xdata = (0:NumConditions-1) * barwidth;
        xdata = xdata - median(xdata) + center(i);
        PSTH(i).xx = (0:0.01:NumConditions-1);
        PSTH(i).xx = PSTH(i).xx - median(PSTH(i).xx) + center(i);
        PSTH(i).yy = zeros(numel(PSTH(i).xx),NumElec);
        options = fitoptions(g);
    else
        
    end
    
    idx = time >= 0 & time < 0.1;
    PSTH(i).stimSC = squeeze(sum(PSTH(i).PSTH(idx,:,:,:),1));
    PSTH(i).stimRate = PSTH(i).stimSC / 0.1;
    idx = time >= -0.2 & time < 0;
    PSTH(i).blankSC = squeeze(sum(PSTH(i).PSTH(idx,:,:,:),1));
    PSTH(i).blankRate = PSTH(i).blankSC / 0.2;
    PSTH(i).meanSC = squeeze(mean(PSTH(i).stimSC,1));
    PSTH(i).meanRate = squeeze(mean(PSTH(i).stimRate,1));
    
    % simple ttest
    PSTH(i).h = ttest(PSTH(i).stimRate,PSTH(i).blankRate,'tail','right');
    idx = isnan(PSTH(i).h);
    PSTH(i).h(idx) = 0;
    PSTH(i).h = squeeze(PSTH(i).h);
    % fit Guassian to x and y to determine the RF location
    if i < 3
        for thisElec = 1:NumElec
            if sum(PSTH(i).h(:,thisElec) == 1) == 0
                % no significantly larger than 0 conditions
                PSTH(i).gfit{thisElec} = [];
                PSTH(i).gof{thisElec} = [];
            else
                ydata = PSTH(i).meanRate(:,thisElec);
                options.Lower = [0.1,min(xdata),barwidth,0];
                options.Upper = [max(ydata),max(xdata),max(xdata)-min(xdata),max(ydata)/2];
                options.startpoint = [max(ydata)/10,(max(xdata)-min(xdata))/2,(max(xdata)-min(xdata))/2,min(ydata)];
                [gfit,gof] = fit(xdata',ydata,g,options);
                PSTH(i).yy(:,thisElec) = gfit(PSTH(i).xx);
                PSTH(i).gfit{thisElec} = gfit;
                PSTH(i).gof{thisElec} = gof;
            end
        end
    end
end

% % SC 32

% x

h=figure('position',[10 200 1600 800]);hold on;
title(pathx)
sessionidx = 1;
for thisElec = 1:PSTH(sessionidx).NumElec
    EID = PSTH(sessionidx).electrodes(thisElec);
    subplot(size(arraymap,1),size(arraymap,2),arraymap(EID)),hold on;
    NumConditions = PSTH(sessionidx).NumConditions;
    xdata = (0:NumConditions-1) * barwidth;
    xdata = xdata - median(xdata) + center(1);
    plot(xdata,PSTH(sessionidx).meanRate(:,thisElec));
    plot(PSTH(sessionidx).xx,PSTH(sessionidx).yy(:,thisElec))
    xlim([min(xdata),max(xdata)])
    title(['Elec ',num2str(EID)])
end
saveas(h,[figurepath,slash,'RFx',sessionname,'.png'],'png');

% y
h=figure('position',[10 200 1600 800]);hold on;
title(pathy)
sessionidx = 2;
for thisElec = 1:PSTH(sessionidx).NumElec
    EID = PSTH(sessionidx).electrodes(thisElec);
    subplot(size(arraymap,1),size(arraymap,2),arraymap(EID)),hold on;
    xlim([0 size(PSTH(sessionidx).meanSC,1)])
    NumConditions = PSTH(sessionidx).NumConditions;
    xdata = (0:NumConditions-1) * barwidth;
    xdata = xdata - median(xdata) + center(2);
    plot(xdata,PSTH(sessionidx).meanRate(:,thisElec));
    plot(PSTH(sessionidx).xx,PSTH(sessionidx).yy(:,thisElec))
    xlim([min(xdata),max(xdata)])
    title(['Elec ',num2str(EID)])
end
saveas(h,[figurepath,slash,'RFy',sessionname,'.png'],'png');

% o
% figure
% sessionidx = 3;
% for thisElec = 1:PSTH(sessionidx).NumElec
%     EID = PSTH(sessionidx).electrodes(thisElec);
%     subplot(6,6,arraymap(EID)),hold on;
%     xlim([0 size(PSTH(sessionidx).meanSC,1)])
%     NumConditions = PSTH(sessionidx).NumConditions;
%     xdata = (0:NumConditions-1) * barwidth;
%     xdata = xdata - median(xdata) + center(2);
%     plot(xdata,PSTH(sessionidx).meanRate(:,thisElec));
% %     plot(PSTH(sessionidx).xx,PSTH(sessionidx).yy(:,thisElec))
%     xlim([min(xdata),max(xdata)])
%     title(['Elec ',num2str(EID)])
% end

% plot RF
fh=figure('position',[10 200 1600 1200]);
hold on
plot([-20 20],[0,0],'blue')
plot([0 0],[-20,15],'blue')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
title([{'Gaussion fitting to RF mapped with grating bars'},{['SessionName = ' sessionname]}])
for thisElec = 1:PSTH(1).NumElec
    % gof
    if isempty(PSTH(1).gof{thisElec}) || isempty(PSTH(2).gof{thisElec})
        continue
    end
    h = PSTH(1).gof{thisElec}.rsquare >= rsquarethres & PSTH(2).gof{thisElec}.rsquare >= rsquarethres;
    elecRF.EID(thisElec) = PSTH(1).electrodes(thisElec);
    elecRF.x(thisElec) = PSTH(1).gfit{thisElec}.b;
    elecRF.xwidth(thisElec) = PSTH(1).gfit{thisElec}.c*xSigma - barwidth;
    elecRF.xRsquare(thisElec) = PSTH(1).gof{thisElec}.rsquare;
    elecRF.y(thisElec) = PSTH(2).gfit{thisElec}.b;
    elecRF.ywidth(thisElec) = PSTH(2).gfit{thisElec}.c*xSigma - barwidth;
    elecRF.yRsquare(thisElec) = PSTH(2).gof{thisElec}.rsquare;
    if h
        EID = PSTH(1).electrodes(thisElec);
        RF.x = PSTH(1).gfit{thisElec}.b;
        RF.wx = PSTH(1).gfit{thisElec}.c * xSigma - barwidth;
        RF.y = PSTH(2).gfit{thisElec}.b;
        RF.wy = PSTH(2).gfit{thisElec}.c * xSigma - barwidth;
        text(RF.x+0.1,RF.y+0.1,num2str(EID))
        % center
        plot(RF.x,RF.y,'k.')
        % left
        plot([RF.x-RF.wx, RF.x-RF.wx],[RF.y-RF.wy,RF.y+RF.wy],'black')
        % top
        plot([RF.x-RF.wx,RF.x+RF.wx],[RF.y+RF.wy,RF.y+RF.wy],'black')
        % right
        plot([RF.x+RF.wx, RF.x+RF.wx],[RF.y-RF.wy,RF.y+RF.wy],'black')
        % bottom
        plot([RF.x-RF.wx,RF.x+RF.wx],[RF.y-RF.wy,RF.y-RF.wy],'black')
    end
    
end
saveas(fh,[figurepath,slash,'RFmap',sessionname,'.png'],'png');
FigurePath.x = [figurepath,slash,'RFx',sessionname,'.png'];
FigurePath.y = [figurepath,slash,'RFy',sessionname,'.png'];
FigurePath.RFmap = [figurepath,slash,'RFmap',sessionname,'.png'];
close all