datasource = "MySQLNative";
username = "root";
password = "matlab";

conn = mysql('127.0.0.1','mmsteam','mms_issi$2026')

%%
opts = databaseConnectionOptions("MySQL");
opts.Server = "127.0.0.1";
opts.Port = 899;               % set if nonstandard
opts.User = "mmsteam";
opts.Password = "mms_issi$2026";
conn = connect(opts);

%%
javaaddpath('/Users/cecilianorgren/MATLAB/cn-matlab/mathworks/mysql-connector-j-9.6.0/mysql-connector-j-9.6.0.jar')

dbname   = 'MMS';
username = 'mmsteam';
password = ''; % mms_issi$2026

conn = database( ...
    dbname, ...
    username, ...
    password, ...
    'Vendor', 'MySQL', ...
    'Server', '127.0.0.1', ...
    'PortNumber', 3308);

conn.Message

%%
url = 'jdbc:mysql://127.0.0.1:3308/MMS';
username = 'mmsteam';
password = '';

conn = java.sql.DriverManager.getConnection(url, username, password);