@ECHO OFF
setlocal
CALL activate magine_36
set PYTHONPATH=%PYTHONPATH%;E:\PycharmProjects\PycharmProjects\Magine
make html
pause
endlocal