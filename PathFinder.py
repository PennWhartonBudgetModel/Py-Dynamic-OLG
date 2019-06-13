##
# Path finder for dynamic model source code, inputs, and outputs.
#
#

import os as path
import warning
import pandas as pd
import datetime
import subprocess
import re

class PathFinder:
    
    scenario
    
    # List of Major versions, by component interface, which the current code will accept
    #     TEMP: During transition to this system, an empty ([]) version
    #     means not versioned.
    majorVersions = {
            'Microsim': {'microsim': [], 'entry_exit_rates': [], 'projections' : []},
            'TaxCalculator': {'taxcalculator': '1'},
            'OASIcalculator': {'oasicalculator': []},
            'DynamicModel': {'calibration': [], 'outofmodel': [], 'openness': '1'}
            }
    
    # Specify inputsets of interface versions, organized by component
    #   Production set is for Production & Development modes
    #    Testing set is for Testing mode
    inputversions = {
            'Production': {
                    'Microsim': {'microsim': '2019-03-21-03-24-jricco-GaleRuns',
                          'entry_exit_rates': '2019-01-31-07-19-njanetos-0c2b55d',
                          'projections': '2019-03-15-00-00-jricco-LatestProjectionsBaseline',
                          'openness': '2019-01-08-09-11-njanetos-transition'},
                    'TaxCalculator': {'taxcalculator': '2018-03-13-00-00-00-LatestBaseline'},
                    'OASIcalculator': {'oasicalculator': '2019-03-18-13-20-ses-2ab2306'},
                    'DynamicModel': {'calibration': '2017-10-20-05-56-efraim-f0a4c45',
                                  'outofmodel': '2019-05-03-13-49-jhuntley-7524391'}},
          'Testing': {
                  'Microsim': {'microsim': '2019-01-08-09-11-njanetos-transition',
                       'entry_exit_rates': '2019-01-31-07-19-njanetos-0c2b55d',
                       'projections': '2019-03-15-00-00-jricco-LatestProjectionsBaseline',
                       'openness': '2019-01-08-09-11-njanetos-transition'},
               'TaxCalculator': {'taxcalculator': '2018-03-13-00-00-00-LatestBaseline'},
               'OASIcalculator': {'oasicalculator': '2019-03-18-13-20-ses-2ab2306'},
               'DynamicModel': {'calibration': '2017-10-20-05-56-efraim-f0a4c45',
                       'outofmodel': '2019-05-03-13-49-jhuntley-7524391'}}}
    
    # Constructor
    def __init__(scenario):
        self.scenario = scenario
        return self
    
    # Get path for Openness interface
    def getOpennessInputPath(self, inputName):
        path = self.getInputFilePath('DynamicModel', 'openness', self.scenario.VersionOpenness, inputName)
        return path
    
    # Get path for OASIcalculator interface
    def getOASIcalculatorInputPath(self, inputName):
        path = self.getInputFilePath('OASIcalculator', 'oasicalculator', self.scenario.VersionOASIcalculator, inputName)
        return path
    
    # Get path for TaxCalculator interface
    def getTaxCalculatorInputPath(self, inputName):
        path = self.getInputFilePath('TaxCalculator', 'taxcalculator', self.scenario.VersionTaxCalculator, inputName)
        return path
    
    # Get path for Microsim interface
    def getMicrosimInputPath(self, inputName):
        path = self.getInputFilePath('Microsim', 'microsim', self.scenario.VersionMicrosim, inputName)
        return path
    
    # Get path for Demographics interface
    def getDemographicsInputPath(self, inputName):
        path = self.getInputFilePath('Microsim', 'entry_exit_rates', self.scenario.VersionEntryExitRates, inputName)
        return path
    
    # Get path for Demographics interface
    def getProjectionsInputPath(self, inputName):
        path = self.getInputFilePath('Microsim', 'projections', self.scenario.VersionProjections, inputName)
        
    # Get path for OutOfModel interface
    def getOutOfModelInputPath(self, inputName):
        path = self.getInputFilePath('DynamicModel', 'outofmodel', self.scenario.VersionOutOfModel, inputName)
        
    # Get versioned location of an input source
    #    Inputs:
    #           interfaceDir    : base of interface, expected location of map file
    #           inputName       : the name of the file e.g., 'aggregates', all
    #           are expected to end in 'csv'
    #   Output: 
    #           inputFile : location of input files
    def getInputFilePath(self, component, interface, version, inputName):
        
        # Add major versionID to interface, validation will occur through file system
        majorVersion    = PathFinder.majorVersions[component][interface]
        interface      += majorVersion
        
        # Construct the path to the input directory
        interfaceDir    = os.path.join(PathFinder.getComponentRootDir(component), 
                              'Interfaces', version, interface)
        
        # Load mapping table from file
        mapfile = path.join(interfaceDir, 'map.csv')
        if path.exists(mapfile):
            # TBD: This is for robustness of mis-named files -- should all
            # be map.csv
            mapfile = path.join(interfaceDir, 'Map.csv')
        
        # TEMP: All components should have mapfiles
        if path.exists(mapfile):
            warning( 'No map file for %s' % interfaceDir)
            path = path.join(interfaceDir, inputName + '.csv')
            
        # TBD: There should be an 'ID' column
        map_df = pd.read_csv(mapfile, index_col = 0, header = None).transpose()
        
        # Specify input parameters to use for matching
        # Assumes a 1-to-1 relationship between input parameters and identically named dynamic model scenario properties
        matchparams = ['OpennessPath', 'TaxCode', 'OutlaysPolicy', 'OutOfModel', 'Description']
        
        # Filter input parameters for those present in mapping table
        matchparams = [v for v in matchparams if v in map_df.columns.values]
        
        # Identify matching input scenarios
        f = lambda input_scenario: np.all([self.scenario[param] == input_scenario[param] for param in matchparams])
        match = {k: f(v) for k, v in map_df.to_dict(orient='list').items()}
        
        # Check for singular match
        if sum(match) < 1:
            for param in matchparams:
                print( 'Scenario.%s="%s" ' % (param, self.scenario[param]))
            print( '\n' )
            raise Exception( 'No matches found for scenario in mapfile %s' % mapfile )
        else:
            if sum(match) > 1:
                print( '\nDuplicate ids: ')
                for i in range(len(match)):
                    if match[i]:
                        s = map_df.columns.values[i]
                        print( ' %s ' % str(s))
                print('\n')
                raise Exception( 'More than one input scenario found corresponding to dynamic model scenario. Mapfile %s' % mapfile)
        
        # Extract ID of matching input scenario
        # TBD: Get by column name "ID" not just first column
        id = map.columns.values[match]
        
        # TBD: All ID files should be in their own directories
        #      For robustness, allow old style of location
        path = os.join.path(interfaceDir, id, inputName + '.csv' )
        if not os.path.isfile(path):
            path = os.join.path(interfaceDir, inputName + '_' + id + '.csv' )
            if not os.path.join(path):
                raise Exception('Cannot find input file: %s' % path)

    
    # Singleton execution mode wrapper
    #   Serves as both getter and setter
    #   Required to work around lack of non-constant static variables in Matlab
    def ExecutionMode(newmode):
        
        global mode_
        
        # Set to new execution mode if provided
        if 'newmode' in globals():
            if newmode in ['Development', 'Production', 'Testing']:
                mode_ = newmode
            else:
                raise Exception('Execution mode must be in (''Development'', ''Testing'', ''Production'')')
        else:
            # Initialize to development mode if uninitialized
            if mode_ == None:
                mode_ = 'Development'
        
        mode = mode_
        
        return mode
    
    
    # Singleton wrapper for InputRootDir
    #   Serves as both getter and setter
    #   Required to work around lack of non-constant static variables in Matlab
    def InputRootDir(newdir):
        
        global dir_
        
        # Set to new execution mode if provided
        if 'newdir' in globals():
            dir_ = newdir
        else:
            # Initialize to HPCC root dir uninitialized
            if dir_ == None:
                dir_ = PathFinder.getHpccRootDir()

        thedir = dir_
        
        return thedir
    
    
    
    # Get HPCC root directory
    #   Local path used if executing on the HPCC or on AWS through the HPCC
    #   Server path used if executing elsewhere
    def getHpccRootDir():
        if PathFinder.isHPCCRun():
            d = os.join.path(os.sep, 'home', 'mnt', 'projects')
        else:
            d = os.join.path(os.sep, os.sep, 'hpcc-ppi.wharton.upenn.edu')
            
        hpccrootdir = os.path.join(d, 'ppi')
    
        return hpccrootdir
    
    # Get component root directory
    #   Defaults to dynamic model if no component specified
    def getComponentRootDir(component):
        if not 'component' in 'var' or component == None:
            component = 'DynamicModel'
        componentrootdir = os.path.join(PathFinder.InputRootDir(), component)
        return componentrootdir
    
    # Get working root directory
    def getWorkingRootDir():
        if PathFinder.ExecutionMode() == 'Development':
            workingrootdir = os.join.path(PathFinder.getSourceDir(), 'Working')
        elif PathFinder.ExecutionMode() == 'Testing':
            workingrootdir = os.join.path(PathFinder.getSourceDir(), 'Working')
        elif PathFinder.ExecutionMode() == 'Production':
            workingrootdir = os.join.path(PathFinder.getComponentRootDir(), 'Internal', PathFinder.getCommitTag())
        
        return workingrootdir
    
    # Get output directory
    def getOutputDir(interface):
        if PathFinder.ExecutionMode() == 'Development':
            outputrootdir = PathFinder.getWorkingRootDir()
        elif PathFinder.ExecutionMode() == 'Testing':
            outputrootdir = PathFinder.getWorkingRootDir()
        elif PathFinder.ExecutionMode() == 'Production':
            outputrootdir = os.join.path(PathFinder.getComponentRootDir(), 'Interfaces', PathFinder.getCommitTag())
        
        outputdir = os.join.path(outputrootdir, interface)
    
        return outputdir

    
    # Get unique identifying tag for active Git commit
    def getCommitTag():
        
        # Get commit date (%cd), author email address (%ae), and abbreviated hash (%h)
        result = subprocess.run('git --no-pager log -1 --format=%cd,%ae,,%h --date=iso', stdout=subprocess.PIPE, shell = True)
        commitlog = result.stdout.decode('utf-8')
        
        # Extract commit date and reformat
        commitdate = re.search('.*? .*?(?= )', commitlog).groups(1)
        commitdate = datetime.strptime(commitdate, 'yyyy-mm-dd-HH-MM')
        
        # Extract author username from email address
        commitauthor = re.search('(?<=,).*?(?=@)', commitlog).groups(1)
        
        # Extract abbreviated hash
        commithash = re.search('(?<=,,).*?(?=\n)', commitlog).groups(1)
        
        # Construct commit identifier
        committag = '%s-%s-%s' % (commitdate, commitauthor, commithash)
        
        return committag

    
    # Check if running on HPCC (or AWS)
    def isHPCCRun():
        b = ~isempty(regexp(getenv('HOSTNAME'), '^(hpcc|aws)', 'once'))
        return b
    
    % Set input root to be HPCC
    function [] = setInputRootToHPCC()
        PathFinder.InputRootDir(PathFinder.getHpccRootDir()); % Kludge to set 'static' variable
    end
    
    % Set input root to be Local
    function [] = setInputRootToLocal()
        mode = PathFinder.ExecutionMode();
        if( strcmp(mode, 'Development')) 
            inputrootdir = fullfile(PathFinder.getSourceDir(), 'CachedInput');
            PathFinder.InputRootDir(inputrootdir);  % Kludge to set 'static' variable
        else
            throw(MException(   'PathFinder:SETINPUTROOTDIR'                                        ...
                                ,   'Cannot use local input cache in Testing or Production modes.'  ...
                                ));
        end
    end
    
    % Get execution mode
    function mode = getExecutionMode()
        mode = PathFinder.ExecutionMode();
    end
    
    % Set execution mode to development
    function [] = setToDevelopmentMode()
        PathFinder.ExecutionMode('Development');
    end
    
    % Set execution mode to Testing
    function [] = setToTestingMode()
        PathFinder.ExecutionMode('Testing');
    end
    
    % Set execution mode to production
    function [] = setToProductionMode()
        
        % Check for uncommitted changes
        [~, uncommitted] = system('git status -s');
        
        if (isempty(uncommitted))
            PathFinder.ExecutionMode('Production');
        else
            error( 'Uncommitted changes found. Cannot set to production mode.' );
        end
        
    end
    
    
    % Get list of input set used for this ExecutionMode
    function [versionset] = getInputSet()
        % Production & Development input sets are the same
        if( strcmp( PathFinder.ExecutionMode(), 'Testing' ) )
            inputset = 'Testing';
        else
            inputset = 'Production';
        end
        
        versionset = {};
        for component = fieldnames(PathFinder.inputversions.(inputset))'
            for interface = fieldnames(PathFinder.inputversions.(inputset).(component{1}))'
                version  = PathFinder.inputversions.(inputset).(component{1}).(interface{1});
                versionset{end+1} = {component{1}, interface{1}, version};
            end
        end
    end %getInputSet
    
    
    % Get source code directory
    function [sourcedir] = getSourceDir()
        sourcedir = fileparts(mfilename('fullpath'));
    end
    
    
    % Get cached output directory for a scenario
    function [cachedir] = getCacheDir(scenario)
        cachedir = fullfile(PathFinder.getWorkingRootDir(), ...
            scenario.basedeftag, scenario.counterdeftag, scenario.transitiontag);
    end
    
    
    % Get a temporary working directory for a scenario
    function [tempDir, tempTag] = getTempWorkingDir(scenario, callerTag)
        processID   = feature( 'getPID' );   % NOTE: This uses an undocumented/unsupported feature
        tempTag     = sprintf( '%u+%s', processID, callerTag );
        tempDir     = fullfile(PathFinder.getWorkingRootDir(),                       ...
                        scenario.basedeftag, scenario.counterdeftag, scenario.transitiontag );
        tempDir     = sprintf('%s%s', tempDir, tempTag );
    end

    
    % Get tag path for identifying scenarios
    function [tagPath] = getScenarioPathTag(scenario)
        tagPath = fullfile(          ...
            scenario.basedeftag,     ...
            scenario.counterdeftag,  ...
            scenario.transitiontag   ...
        );
    end
    
    % Get output directory
    function [outputdir] = getCalibrationOutputDir()        , outputdir = PathFinder.getOutputDir('calibration' ); end
    function [outputdir] = getSeriesOutputDir()             , outputdir = PathFinder.getOutputDir('series'      ); end
    function [outputdir] = getTransitionMatrixOutputDir()   , outputdir = PathFinder.getOutputDir('transition'  ); end
    
    
    
end


end
