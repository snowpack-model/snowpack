:: Copy and edit this file for your own simulations!
@ECHO OFF
ECHO Running Snowpack for %USERNAME% on %DATE% %TIME%

snowpack -c cfgfiles/io_res1exp.ini -e 1996-06-17T00:00

:: wait for the user to press a key before closing the window
pause
