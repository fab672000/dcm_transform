@echo off
REM
REM dcm_transform example of use below:
REM in order to use that example you need to:
REM 	1. add any DICOM image file tree structure under the 'original' sub-directory where this batch executes
REM		2. run this script under the windows command line prompt or double-click on it from a Windows explorer window.
REM

set rootDir=%~dp0

set inputDir=original
set outputDir=transformed

set script=dcm_transform.py
set options=-r -an anon

pushd .

cd %rootDir%

REM run recursively through any dicom files contained in the inputDir tree and generate an ouputTree with the resulting edited dicom files
python %script% %options% -sn 66 -x +10 -y +20 -z -40 %inputDir% %outputDir%

rem end-up with a new regeneration of original with same anon id:
python %script% %options% -sn 64 %inputDir% new-original 

popd

pause