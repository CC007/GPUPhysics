@echo off
setlocal ENABLEDELAYEDEXPANSION
set backPath=
if "%4"=="" goto folderNotSet
if not "%4"=="" goto folderSet

:folderNotSet
set /P folder=Enter path to source files: %=%
goto slashTest

:folderSet
set folder=%4
goto slashTest 


:slashTest
if not "%folder:~-1%"=="/" goto addFolderSlash
goto countSlashes

:addFolderSlash
set folder=%folder%/
goto countSlashes

:countSlashes
call charCount %folder% /
set slashCount=%errorlevel%
for /l %%x in (1,1,%slashCount%) do set backPath=../!backPath!
goto checkSources

:checkSources
if not "%5"=="" goto sourcesSet
set cuSources=main.cu
goto continue

:sourcesSet
set cuSources=%5
goto continue

:continue
echo echo cd %folder% ^> make.txt > make.txt
echo echo make CUSOURCES=\'%cuSources%\' all -f %backPath%Makefile.dat ^>^> make.txt >> make.txt
echo putty-0.63/putty-0.63/unix/plink %1@srv016 -pw %3 -m make.txt >> make.txt
plink %1@linux-login.kvi.nl -pw %2 -m make.txt
pause
