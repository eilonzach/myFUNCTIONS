mail = 'eilonzach@gmail.com';
password = 'grtm3RMI';
server = 'smtp.gmail.com';
props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.port','587');
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Apply prefs and props
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server',server);
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

% Send the email
sendmail('eilonzach@gmail.com', 'Test', 'Msg from MATLAB');
