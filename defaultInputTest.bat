@echo off
setlocal
set "Input="
set "Default=main.cu"
title Input
CScript //E:JScript //Nologo "DefaultInput.bat" "Input" "%Default%"
set /p "Input=> Prompt: "
if defined Input set "Input=%Input:"=%"
echo(%Input%
endlocal
exit /b 0