Python model based on Matlab version of June 24, 2019

In order to run the code, the letter for the drive where HPCC is mounted on the local workstation needs to be hardcoded in pathFinderModule function getHpccRootDir().
Matlab could figure out the drive itself from the HPCC link but I didn't find a way to replicate that in Python.

Code can be run with test parameters by executing the file testSolverModule.

initialGuessDataModule can be used to generate a InitialGuess.pkl file with Scenario objects and Market and Dynamic dictionaries that replicates the InitialGuess.mat file from Matlab. 

To do: the code needs 1) to be sped up 2) to be updated to correspond with the latest Matlab commit.

If you have any questions about this code, do not hesitate to contact Alex Zanca at zanca@sas.upenn.edu.
