@echo off
setlocal ENABLEDELAYEDEXPANSION
set theLine=%1
set theChar=%2
call :strLen "!theLine!"
set totalChars=%errorlevel%
set strippedLine=!theLine:%theChar%=!
call :strLen "!strippedLine!"
set /A n=totalChars-%errorlevel%
exit /B %n%

:strLen
echo "%~1"> StrLen
for %%a in (StrLen) do set /A StrLen=%%~Za-4
del StrLen
exit /B %strLen%