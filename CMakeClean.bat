@echo off
REM windows batch

set MODE=%1%
if "%MODE%"=="" (
    set MODE=quick
)

if "%MODE%"=="quick" (
	for /d /r . %%d in (CMakeFiles) do (
		if exist "%%d" rd /s/q "%%d"
	)
	del /s /q /f CMakeCache.txt
	del /s /q /f cmake_install.cmake
	del /s /q /f Makefile
	del /s /q /f *~
	del /s /q /f *.obj	
)

if "%MODE%"=="deep" (
	for /d /r . %%d in (CMakeFiles) do (
		if exist "%%d" rd /s/q "%%d"
	)
	del /s /q /f CMakeCache.txt
	del /s /q /f cmake_install.cmake
	del /s /q /f Makefile
	del /s /q /f *.obj
	del /s /q /f *.a
	del /s /q /f *.dll
	del /s /q /f *.dll.*
	del /s /q /f *.mod
	for /d /r . %%d in (bin\*) do (
		if exist "%%d" rd /s/q "%%d"
	)
)
