@echo off
REM Compiles the RCS sizing tool with MSVC (Visual Studio 2022) into build\rcs_sizing.exe
setlocal
cd /d "%~dp0"
if not exist build mkdir build
if not exist results mkdir results

where cl >nul 2>nul
if errorlevel 1 call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul

set SRC=src\main.cpp src\layoutAnalysis.cpp src\dvBudget.cpp src\proximityChecks.cpp src\thrustWindow.cpp src\geometrySweeps.cpp src\comparison.cpp
cl /nologo /std:c++17 /EHsc /O2 /W3 /Iinclude /Fobuild\ /Fe:build\rcs_sizing.exe %SRC% || exit /b 1
echo Build OK