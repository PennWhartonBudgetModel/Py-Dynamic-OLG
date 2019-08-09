##
# Path finder for dynamic model source code, inputs, and outputs.
#
#
from scenarioModule import Scenario

import os
import warnings
import pandas as pd
import numpy as np
import datetime
import subprocess
import re
import socket

class PathFinder:
    
    scenario = None
    
    # List of Major versions, by component interface, which the current code will accept
    #     TEMP: During transition to this system, an empty ([]) version
    #     means not versioned.
    majorVersions = {
            'Microsim': {'microsim': r'', 'demographics': r'', 'projections' : r''},
            'TaxCalculator': {'taxcalculator': r'1'},
            'OASIcalculator': {'oasicalculator': r'1'},
            'DynamicModel': {'calibration': r'', 'outofmodel': r'', 'laborshock': r'3', 'openness': r'1'}
            }
    
    # Specify inputsets of interface versions, organized by component
    #   Production set is for Production & Development modes
    #    Testing set is for Testing mode
    inputversions = {
            'Production'    : {
                    'Microsim'  : {'microsim'         : '2019-03-21-03-24-jricco-GaleRuns',                  
                                   'demographics'     : '2019-05-10-15-00-jhuntley-6f1d22',                  
                                   'projections'      : '2019-03-15-00-00-jricco-LatestProjectionsBaseline',
                                   'openness'         : '2019-01-08-09-11-njanetos-transition' },               
                'TaxCalculator' : {'taxcalculator'    : '2018-03-13-00-00-00-LatestBaseline'},           
               'OASIcalculator' : {'oasicalculator'   : '2019-03-18-13-20-ses-2ab2306'},           
                  'DynamicModel': {'calibration'      : '2017-10-20-05-56-efraim-f0a4c45',
                                   'outofmodel'       : '2019-05-03-13-49-jhuntley-7524391',                 
                                   'laborshock'       : '2019-05-03-13-49-jhuntley-7524391'}},               
            'Testing'       : {
                 'Microsim'      : {'microsim'         : '2019-01-08-09-11-njanetos-transition',             
                                    'demographics'     : '2019-05-10-15-00-jhuntley-6f1d22',                  
                                    'projections'      : '2019-03-15-00-00-jricco-LatestProjectionsBaseline', 
                                    'openness'         : '2019-01-08-09-11-njanetos-transition' },               
                 'TaxCalculator' : {'taxcalculator'    : '2018-03-13-00-00-00-LatestBaseline'   },           
                 'OASIcalculator': {'oasicalculator'   : '2019-03-18-13-20-ses-2ab2306'         },          
                 'DynamicModel'  : {'calibration'      : '2017-10-20-05-56-efraim-f0a4c45',           
                                    'outofmodel'       : '2019-05-03-13-49-jhuntley-7524391',                
                                    'laborshock'       : '2019-05-03-13-49-jhuntley-7524391',    }}}           
    
    # Constructor
    def __init__(self, scenario):

        self.scenario = scenario

    
    # Get path for Openness interface
    def getOpennessInputPath(self, inputName):
        path = self.getInputFilePath('DynamicModel', r'openness', self.scenario.VersionOpenness, inputName)
        return path
    
    # Get path for OASIcalculator interface
    def getOASIcalculatorInputPath(self, inputName):
        path = self.getInputFilePath('OASIcalculator', r'oasicalculator', self.scenario.VersionOASIcalculator, inputName)
        return path
    
    # Get path for TaxCalculator interface
    def getTaxCalculatorInputPath(self, inputName):
        path = self.getInputFilePath('TaxCalculator', r'taxcalculator', self.scenario.VersionTaxCalculator, inputName)
        return path
    
    # Get path for Microsim interface
    def getMicrosimInputPath(self, inputName):
        path = self.getInputFilePath('Microsim', r'microsim', self.scenario.VersionMicrosim, inputName)
        return path
    
    # Get path for Demographics interface
    def getDemographicsInputPath(self, inputName):
        path = self.getInputFilePath('Microsim', r'demographics', self.scenario.VersionDemographics, inputName)
        return path
    
    # Get path for Demographics interface
    def getProjectionsInputPath(self, inputName):
        path = self.getInputFilePath('Microsim', r'projections', self.scenario.VersionProjections, inputName)
        return path
        
    # Get path for OutOfModel interface
    def getOutOfModelInputPath(self, inputName):
        path = self.getInputFilePath('DynamicModel', r'outofmodel', self.scenario.VersionOutOfModel, inputName)
        return path
    
    # Get path for LaborShock interface
    def getLaborShockInputPath(self, inputName):
        path = self.getInputFilePath('DynamicModel', r'laborshock', self.scenario.VersionLaborShock, inputName)
        return path
    
    # Get path for named output location
    def getNamedOutputPath(self, outputName):
        path = os.path.join(PathFinder.getOutputDir( 'named' ), outputName )
        return path
        
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
        interface += majorVersion
        
        # Construct the path to the input directory
        interfaceDir    = os.path.join(PathFinder.getComponentRootDir(component), 
                              r'Interfaces', version, interface)
        
        # Load mapping table from file
        mapfile = os.path.join(interfaceDir, r'map.csv')
        if not os.path.exists(mapfile):
            # TBD: This is for robustness of mis-named files -- should all
            # be map.csv
            mapfile = os.path.join(interfaceDir, r'Map.csv')
            if not os.path.exists(mapfile):
                raise Exception ('Cannot find map file: %s' % mapfile)
        
        '''
        # TEMP: All components should have mapfiles
        if not os.path.exists(mapfile):
            warnings.warn( 'No map file for %s' % interfaceDir)
            path = os.path.join(interfaceDir, inputName + '.csv')
            return path
            '''
            
        # TBD: There should be an 'ID' column
        map_df = pd.read_csv(mapfile, index_col = 0, header = 0)
        
        # Specify input parameters to use for matching
        # Assumes a 1-to-1 relationship between input parameters and identically named dynamic model scenario properties
        scenarioParams = self.scenario.getMapFileMatchParams();
        
        # Filter input parameters for those present in mapping table
        matchparams = scenarioParams.keys()
        matchparams = [v for v in matchparams if v in map_df.columns.values]
        
        # Identify matching input scenarios
        f = lambda input_scenario: np.all([scenarioParams[param] == input_scenario[param] for param in matchparams])
        match = [f(v) for v in map_df.to_dict(orient='records')]
        
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
                        s = map_df.index.values[i]
                        print( ' %s ' % str(s))
                print('\n')
                raise Exception( 'More than one input scenario found corresponding to dynamic model scenario. Mapfile %s' % mapfile)
        
        # Extract ID of matching input scenario
        # TBD: Get by column name "ID" not just first column
        id = map_df.index.values[match]
        id = r''.join(id)
        
        # TBD: All ID files should be in their own directories
        #      For robustness, allow old style of location
        path = os.path.join(interfaceDir, id, inputName + r'.csv' )
        if not os.path.isfile(path):
            path = os.path.join(interfaceDir, inputName + r'_' + id + r'.csv' )
            if not os.path.join(path):
                raise Exception('Cannot find input file: %s' % path)
                
        return path

    
    # Singleton execution mode wrapper
    #   Serves as both getter and setter
    #   Required to work around lack of non-constant static variables in Matlab
    # TBD find something better than mutable default arguments to preserve mode value beteween static calls
    @staticmethod
    def ExecutionMode(newmode = None, mode = []):
        
        if newmode != None:
            if newmode in ['Development', 'Production', 'Testing']:
                if len(mode) != 0:
                    mode.pop()
                mode.append(newmode)
            else:
                raise Exception('Execution mode must be in (''Development'', ''Testing'', ''Production'')')
        elif len(mode) == 0:
            mode.append('Development')

        return mode[0]
    
    # Singleton wrapper for InputRootDir
    #   Serves as both getter and setter
    #   Required to work around lack of non-constant static variables in Matlab
    # TBD find something better than mutable default arguments to preserve mode value beteween static calls
    @staticmethod
    def InputRootDir(newdir = None, dir = []):
        
        if newdir != None:
            if len(dir) != 0:
                dir.pop()
            dir.append(newdir)
        
        elif len(dir) == 0:
            dir.append(PathFinder.getHpccRootDir())
        
        return dir[0]
    
    # Get HPCC root directory
    #   Local path used if executing on the HPCC or on AWS through the HPCC
    #   Server path used if executing elsewhere
    @staticmethod
    def getHpccRootDir():
        if PathFinder.isHPCCRun():
            d = os.path.join(os.sep, 'home', 'mnt', 'projects')
        else:
            #TBD won't work if hpcc is not mounted as Z:
            d = r'Z:'
    
        return d
    
    # Get component root directory
    #   Defaults to dynamic model if no component specified
    @staticmethod
    def getComponentRootDir(component=r'DynamicModel'):
        componentrootdir = os.path.join(PathFinder.InputRootDir(), component)
        return componentrootdir
    
    # Get working root directory
    @staticmethod
    def getWorkingRootDir():
        if PathFinder.ExecutionMode() == 'Development':
            workingrootdir = os.path.join(PathFinder.getSourceDir(), 'Working')
        elif PathFinder.ExecutionMode() == 'Testing':
            workingrootdir = os.path.join(PathFinder.getSourceDir(), 'Working')
        elif PathFinder.ExecutionMode() == 'Production':
            workingrootdir = os.path.join(PathFinder.getComponentRootDir(), 'Internal', PathFinder.getCommitTag())
        
        return workingrootdir
    
    # Get output directory
    @staticmethod
    def getOutputDir(interface):
        if PathFinder.ExecutionMode() == 'Development':
            outputrootdir = PathFinder.getWorkingRootDir()
        elif PathFinder.ExecutionMode() == 'Testing':
            outputrootdir = PathFinder.getWorkingRootDir()
        elif PathFinder.ExecutionMode() == 'Production':
            outputrootdir = os.path.join(PathFinder.getComponentRootDir(), 'Interfaces', PathFinder.getCommitTag())
        
        outputdir = os.path.join(outputrootdir, interface)
    
        return outputdir

    
    # Get unique identifying tag for active Git commit
    @staticmethod
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
    @staticmethod
    def isHPCCRun():
        host = socket.gethostname()
        b = re.search('^(hpcc|aws)', host)
        return (b != None)
    
    # Set input root to be HPCC
    @staticmethod
    def setInputRootToHPCC():
        PathFinder.InputRootDir(PathFinder.getHpccRootDir()) # Kludge to set 'static' variable
    
    # Set input root to be Local
    @staticmethod
    def setInputRootToLocal():
        mode = PathFinder.ExecutionMode()
        if mode == 'Development': 
            inputrootdir = os.path.join(PathFinder.getSourceDir(), 'CachedInput')
            PathFinder.InputRootDir(inputrootdir)  # Kludge to set 'static' variable
        else:
            raise Exception('Cannot use local input cache in Testing or Production modes.')
    
    # Get execution mode
    @staticmethod
    def getExecutionMode():
        mode = PathFinder.ExecutionMode()
        return mode
    
    # Set execution mode to development
    @staticmethod
    def setToDevelopmentMode():
        PathFinder.ExecutionMode('Development')
    
    # Set execution mode to Testing
    @staticmethod
    def setToTestingMode():
        PathFinder.ExecutionMode('Testing')
    
    # Set execution mode to production
    @staticmethod
    def setToProductionMode():
        
        # Check for uncommitted changes
        result = subprocess.run('git status -s', stdout=subprocess.PIPE, shell = True)
        uncommitted = result.stdout.decode('utf-8')
        
        if len(uncommitted) == 0:
            PathFinder.ExecutionMode('Production')
        else:
            raise Exception( 'Uncommitted changes found. Cannot set to production mode.' )
    
    
    # Get list of input set used for this ExecutionMode
    @staticmethod
    def getInputSet():
        # Production & Development input sets are the same
        if PathFinder.ExecutionMode() == 'Testing' :
            inputset = 'Testing'
        else:
            inputset = 'Production'
        
        versionset = np.zeros((0,3))
        for component in PathFinder.inputversions[inputset].keys():
            for interface in PathFinder.inputversions[inputset][component].keys():
                version = PathFinder.inputversions[inputset][component][interface]
                versionset = np.vstack((versionset, np.array([component, interface, version])))

        return versionset
    
    # Get source code directory
    @staticmethod
    def getSourceDir():
        fullpath = os.path.realpath(__file__)
        (drive, path) = os.path.splitdrive(fullpath)
        (path, file)  = os.path.split(path)
        # return (drive, path, file)
        return os.path.join(drive, path)
    
    # Get cached output directory for a scenario
    @staticmethod
    def getCacheDir(scenario):
        cachedir = os.path.join(PathFinder.getWorkingRootDir(),
            scenario.basedeftag, scenario.counterdeftag, scenario.transitiontag)
        return cachedir
    
    
    # Get a temporary working directory for a scenario
    @staticmethod
    def getTempWorkingDir(scenario, callerTag):
        processID   = os.getpid()   # different from Matlab
        tempTag     = '%u+%s' % (processID, callerTag)
        tempDir     = os.path.join(PathFinder.getWorkingRootDir(),                       
                        scenario.basedeftag, scenario.counterdeftag, scenario.transitiontag )
        tempDir     = '%s%s' % (tempDir, tempTag)
        return (tempDir, tempTag)


    # Get tag path for identifying scenarios
    @staticmethod
    def getScenarioPathTag(scenario):
        tagPath = os.path.join(          
            scenario.basedeftag,     
            scenario.counterdeftag,
            scenario.transitiontag)
        return tagPath
    
    # Get output directory
    @staticmethod
    def getCalibrationOutputDir():
        outputdir = PathFinder.getOutputDir('calibration')
        return outputdir
    @staticmethod
    def getSeriesOutputDir():
        outputdir = PathFinder.getOutputDir('series')
        return outputdir
    @staticmethod
    def getTransitionMatrixOutputDir():
        outputdir = PathFinder.getOutputDir('transition')
        return outputdir
    
