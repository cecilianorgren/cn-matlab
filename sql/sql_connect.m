%% Connect to database, need to do the tunnel thing in terminal first
mysqlcheck = mysql('status');
if (mysqlcheck > 0)
    mysql('open', '127.0.0.1:3308','mmsteam','');
end
clear mysqlcheck;
mysql('use', 'MMS');        