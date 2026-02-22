% In terminal: set up tunnel
% ssh -p899 -L 3308:localhost:3306 mmsteam@81.169.221.160 
% password: mms_issi$2026

% Download JAVA driver
% https://dev.mysql.com/downloads/connector/j/
% I put mine in my matlab directory

% Add to path, I out this is my startup.m file
javaaddpath('/Users/cecilianorgren/MATLAB/cn-matlab/mathworks/mysql-connector-j-9.6.0/mysql-connector-j-9.6.0.jar')


% Enter details for database
dbname   = 'MMS'; % database name
username = 'mmsteam'; % database user name
password = ''; % no password for database

conn = database( ...
    dbname, ...
    username, ...
    password, ...
    'Vendor', 'MySQL', ...
    'Server', '127.0.0.1', ...
    'PortNumber', 3308);

conn.Message % empty if no error
isopen(conn)

%% Some queries
fetch(conn, 'SELECT DATABASE()')
fetch(conn, 'SHOW TABLES')

%%
query = 'SELECT * from GP_tmp';
data = fetch(conn,query);

%%
query = "SELECT X,Y,Z from GP_tmp WHERE shear>45 AND Bup>20 AND FlagStr REGEXP '(^|,)mp(,|$)';";
data = fetch(conn,query);

%%
data_cell = cellfun(@(x)str2num(x),data.shear,'UniformOutput',false);
data_vec = cat(1,data_cell{:});

histogram(data_vec,45:2:180)

%% Find magnetopause events

flag_mp = find_flag(data.FlagStr,'mp');

% Shear for those events
shear_cell = cellfun(@(x) str2num(x),data.shear,'UniformOutput',false);
shear_vec = cat(1,shear_cell{:});

% Beta_up for those events
bup_cell = cellfun(@(x) str2num(x),data.beta_up,'UniformOutput',false);
%bup_vec = cat(1,bup_cell{:});
bup_vec = [bup_cell{:}]';

idx_mp = find(flag_mp);
idx_shear45 = find(shear_vec>45);
idx_bup = find(bup_vec>20);

idx = intersect(idx_mp,idx_shear45);
idx = intersect(idx,idx_bup);

h = setup_subplots(1,2);
isub = 1;

hca = h(isub); isub = isub + 1;
data_tmp = bup_vec(idx);
histogram(hca,data_tmp,0:.5:20)

hca = h(isub); isub = isub + 1;
data_tmp = shear_vec(idx);
histogram(hca,data_tmp,45:5:180)

%%
query = "SELECT X,Y,Z from GP_tmp WHERE shear>45 AND Bup>20 AND FlagStr REGEXP '(^|,)mp(,|$)';";
data = fetch(conn,query);

X = cellfun(@(x)str2num(x),data.X,'UniformOutput',false);
X = cat(1,X{:});

Y = cellfun(@(x)str2num(x),data.Y,'UniformOutput',false);
Y = cat(1,Y{:});

Z = cellfun(@(x)str2num(x),data.Z,'UniformOutput',false);
Z = cat(1,Z{:});

x_edges = 0:1:20;
y_edges =-15:1:15;
[N] = histcn([X, Y],x_edges,y_edges);

localuser = '';
data_R = load(['/Users/' localuser '/Data/MMS/DB_Lalti/Proba_full_mms1_pos.mat']);


