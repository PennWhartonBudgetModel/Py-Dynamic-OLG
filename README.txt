Python model based on Matlab version of June 24, 2019

Parts that have been checked so far and should behave like Matlab:
- convergenceModule
- dogeCobtrolledModule
- dynamicFirmGEModule
- firmModule
- socialSecurityModule
- paramGeneratorModule
- firmTesterModule
- initialGuessModule
- pathFinderModule
- govtDebtModule
- scenarioModule (but still need to check the writeTransitionMatrix function)

modelSolverModule is checked only up to the first call of generate_aggregates and in generate_aggregates up to line 354 "q = lambda F: np.sum(..."

All other files are still being debugged

TBD: replace all references to .mat files (and how these are loaded and saved) with something Python specific, probably .pkl. The scipy library for processing .mat files seems buggy.

In order to run the code, the letter for the drive where HPCC is mounted on the local workstation needs to be hardcoded in pathFinderModule function getHpccRootDir().
Matlab could figure out the drive itself from the HPCC link but I didn't find a way to replicate that in Python (yet).

Code can be run with test parameters by executing the file testSolverModule.

initialGuessDataModule can be used to generate a InitialGuess.pkl file with Scenario objects and Market and Dynamic dictionaries that replicates the InitialGuess.mat file from Matlab. 