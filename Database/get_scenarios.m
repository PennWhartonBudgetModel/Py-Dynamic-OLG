function scenarios = get_scenarios()

% Add JDBC driver to Matlab Java path
javaaddpath('sqljdbc41.jar');

% Establish connection with development database
connection = database('second_chart', 'pwbm', 'HbXk86rabjehD2AN', ...
                      'Vendor', 'Microsoft SQL Server', 'AuthType', 'Server', ...
                      'Server', 'ppi-slcsql.wharton.upenn.edu', 'PortNumber', 49170);

% Get scenarios from Scenario table
o = exec(connection, 'SELECT * FROM Scenario');
scenarios = o.fetch().Data;
close(o);

% Close database connection
close(connection);

end