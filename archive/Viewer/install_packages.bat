REM https://www.howtogeek.com/266621/how-to-make-windows-10-accept-file-paths-over-260-characters/

set root=%userprofile%

call %root%\anaconda3\Scripts\activate.bat base

REM call conda create --name viewer_env python=3.7
call conda env create -f viewer_env.yml python=3.7

call conda activate viewer_env

echo Installation finished successfully.
pause
