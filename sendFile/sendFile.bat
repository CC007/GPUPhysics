@echo off
title Send a folder to the GPU server
if "%4"=="" goto folderNotSet
if not "%4"=="" goto folderSet

:folderNotSet
set defaultFolder=files
set defaultFolderPath=../projects
CScript //E:JScript //Nologo "DefaultInput.bat" "Send a folder to the GPU server" "%defaultFolder%"
set /P folder=Enter folder name: %=%
goto pathNotSet

:folderSet
set folder=%4
if not "%5"=="" goto pathSet
goto pathNotSet

:pathNotSet
CScript //E:JScript //Nologo "DefaultInput.bat" "Send a folder to the GPU server" "%defaultFolderPath%"
set /P folderPath=Enter path to folder (leave empty if the folder is located at the current location): %=%
if not "%folderPath%"=="" goto slashTest
goto continue

:pathSet
set folderPath=%5
goto slashTest

:slashTest
if not "%folderPath:~-1%"=="\" goto addFolderPathSlash
goto continue

:addFolderPathSlash
set folderPath=%folderPath%\
goto continue


:continue
echo rm -r %folder% > delFiles.txt
echo put -r %folderPath%%folder% > winToKvi.txt
plink %1@linux-login.kvi.nl -pw %2 -m delFiles.txt
psftp %1@linux-login.kvi.nl -pw %2 -b winToKvi.txt
echo echo rm -r %folder%^>delFiles.txt > kvicommands.txt
echo echo put -r %folder%^>kviToServer.txt >> kvicommands.txt
echo putty-0.63/putty-0.63/unix/plink %1@srv016 -pw %3 -m delFiles.txt >> kvicommands.txt
echo putty-0.63/putty-0.63/unix/psftp %1@srv016 -pw %3 -b kviToServer.txt >> kvicommands.txt
plink %1@linux-login.kvi.nl -pw %2 -m kvicommands.txt
pause