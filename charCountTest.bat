@echo off
setlocal ENABLEDELAYEDEXPANSION
call charCount.bat blablabla b
set slashCount=%errorlevel%
set backPath=
for /l %%x in (1,1,%slashCount%) do set backPath=../!backPath!
)
echo %backPath%
pause