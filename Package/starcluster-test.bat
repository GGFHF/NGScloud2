@echo off

rem ----------------------------------------------------------------------------

rem This software has been developed by:
rem
rem     GI Sistemas Naturales e Historia Forestal (formerly known as GI Genetica, Fisiologia e Historia Forestal)
rem     Dpto. Sistemas y Recursos Naturales
rem     ETSI Montes, Forestal y del Medio Natural
rem     Universidad Politecnica de Madrid
rem     https://github.com/ggfhf/
rem
rem Licence: GNU General Public Licence Version 3.

rem ----------------------------------------------------------------------------

rem This script executes actions on a cluster in the test environment.

rem ----------------------------------------------------------------------------

rem Set run environment

setlocal EnableDelayedExpansion

set ERROR=0

set PYTHONPATH=.

set STARCLUSTER_CONFIG=.\config\test-ngscloud2-config.txt

rem ----------------------------------------------------------------------------

rem Execute the program starcluster.exe

:STARCLUSTER

starcluster.exe %*
if %ERRORLEVEL% neq 0 (set RC=%ERRORLEVEL% & set ERROR=1 & goto END)

rem ----------------------------------------------------------------------------

:END

if %ERROR% equ 0 (
    rem exit 0
)

if %ERROR% equ 1 (
    echo *** ERROR: StarCluster ended with return code %RC%.
    rem exit %RC%
)

rem ----------------------------------------------------------------------------
