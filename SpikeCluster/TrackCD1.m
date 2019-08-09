function Status = TrackCD(Track)
% Give feedback via email at 10% chance
% Tell the developer who are using the package

Status = 0;
if Track~=1
    return;
end
try
    SCAddress = 'spikeclustershared@gmail.com';
    SCPassword = 'sc7389131';
    RawInternet = getpref('Internet');
    setpref('Internet', 'E_mail', SCAddress);
    setpref('Internet', 'SMTP_Server', 'smtp.gmail.com');
    setpref('Internet', 'SMTP_Username', SCAddress);
    setpref('Internet', 'SMTP_Password', SCPassword);
    Props = java.lang.System.getProperties;
    RawProps{1} = Props.getProperty('mail.smtp.auth');
    RawProps{2} = Props.getProperty('mail.smtp.socketFactory.class');
    RawProps{3} = Props.getProperty('mail.smtp.socketFactory.port');
    Props.setProperty('mail.smtp.auth','true');
    Props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    Props.setProperty('mail.smtp.socketFactory.port','65');
    try
        ComputerName = char(java.net.InetAddress.getLocalHost);
    catch
        [Status, ComputerName] = system('hostname');
        if Status~=0
            if ispc
                ComputerName = getenv('COMPUTERNAME');
            else
                ComputerName = getenv('HOSTNAME');
            end
        end
    end
    ComputerName = strtrim(ComputerName);
    sendmail('spikecluster@gmail.com', ComputerName, '');
    % Restore default settings
    rmpref('Internet');
    if ~isempty(RawInternet)
        FieldName = fieldnames(RawInternet);
        for i = 1:numel(FieldName)
            setpref('Internet', FieldName{i}, RawInternet.(FieldName{i}));
        end
    end
    Props.setProperty('mail.smtp.auth', char(RawProps{1}));
    Props.setProperty('mail.smtp.socketFactory.class', char(RawProps{2}));
    Props.setProperty('mail.smtp.socketFactory.port', char(RawProps{3}));
end