function scenarios = get_scenarios()

% Add JDBC driver to Matlab Java path
javaaddpath('sqljdbc41.jar');

% Establish connection with development database
connection = database('second_chart', 'pwbm', 'HbXk86rabjehD2AN', ...
                      'Vendor', 'Microsoft SQL Server', 'AuthType', 'Server', ...
                      'Server', 'ppi-slcsql.wharton.upenn.edu', 'PortNumber', 49170);

% Get scenarios from Scenario table
o = connection.exec('SELECT * FROM Scenario');
scenarios = o.fetch().Data;
o.close();

% Close database connection
connection.close();

end