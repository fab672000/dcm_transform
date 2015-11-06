@echo off
REM
REM dcm_transform Region of Interest's example of use below:
REM

set rootDir=%~dp0

REM simple squared ROI
python dcm_transform.py single.dcm result.dcm -roi 62 62 4 1023


REM draws 2 rectangles with the same command (6 parameters per rect)
python dcm_transform.py result.dcm result.dcm -rect  82 32 30 20 3 1023 1.0 22 22 20 40 3 1023 1.0

REM draws  ellipse
python dcm_transform.py result.dcm result.dcm -elp 63.5 63.5 20 20 1023 1.0 1 

REM draws 2 filled  rectangles
python dcm_transform.py result.dcm result.dcm -frect 50 10 20 40 500 0.5 80 80 30 30 500 .3

REM draws 2 crosshairs
python dcm_transform.py result.dcm result.dcm -crosshair 30 82 5 1 1024 1.0 60 92 10 2 1024 1.0

REM draws pixels 
python dcm_transform.py result.dcm result.dcm -pixel 30 12 1024 1.0 60 92 500 0.5

pause
start result.dcm
