@ECHO OFF
setlocal
CALL activate magine_37
set PYTHONPATH=%PYTHONPATH%;C:\PycharmProjects\Magine
make html
pause
endlocal