function WPartifactRemoval(instanceinfo,WPN)

SampleRate = instanceinfo(1).samplerate;

% build filters
Fl=1500;
Fn = SampleRate/2;
N = 2;
[B, A] = butter(N,Fl/Fn,'low'); % high pass 500Hz

% Fl=1500;
% F2= 3000;
% Fn = SampleRate/2;
% N = 4;
% [B, A] = butter(N,[Fl/Fn F2/Fn],'stop'); % high pass 500Hz

if mod(WPN,2) == 0
    WPN = WPN + 1;
end
testlength = 2000; % samples
% load electrode data one channel at a time
NumInstance = numel(instanceinfo);
fitmodel = fittype(@(a,b,c,d,x) a + b * x + c * x.^2 + d * x.^3);
fitOptions = fitoptions(fitmodel);
fitOptions.StartPoint = zeros(4,1);
numruns = 1;
for thisInstance = 1:NumInstance
    fprintf('\n')
    for thisElec = 1:128
        elecfp = fopen(instanceinfo(thisInstance).electrodeCachePath{thisElec},'r');
        tmpdata = fread(elecfp,inf,'int16=>double');
        fclose(elecfp);
        %         % high pass filt the data
        
        % remove artifact
        % first test on a piece of the data
        tmpdata = tmpdata(6050000: (6050000 +testlength - 1));
        testAR = zeros(testlength,1);
        fprintf('*')
        figure
        hold on
        plot(tmpdata)
        for thisRun = 1:numruns
            if thisRun == 1
            for i = 1:testlength
                lb = max(1,i-floor(WPN/2));
                ub = min(i + floor(WPN/2),numel(tmpdata));
                xdata = lb:ub;
                x = xdata - i;
                x = x';
                ydata = tmpdata(xdata);
                [fitobject,gof] = fit(x,ydata,fitmodel,fitOptions);
                testAR(i) = fitobject(0);
                if mod(i,500) == 1
                    fprintf('.')
                end
            end
            tmpdata = tmpdata - testAR;
            else
%             plot(testAR)
            tmpdata = filtfilt(B,A,tmpdata);
            end
            
            plot(tmpdata)
        end
        %
        
    end
end
end


% function y = LocalPolynomial(fd,a,x)
% if numel(a) ~= fd + 1
%     error('number of parameters and degree do not match')
% end
% cmd = 'y =';
% for i = 1:(fd + 1)
%     if i == 1
%     else
%         cmd = [cmd, ' + '];
%     end
%     cmd = [cmd, 'a(', num2str(i), ') * x ^ ', num2str(i-1)];
% end
% cmd = [cmd, ';'];
% eval(cmd);
% end
% 
% function err = SumSquareError(Y,y)
%     err = sum((Y - y).^2);
% end