set root=%userprofile%

call %root%\anaconda3\Scripts\activate.bat viewer_env

call conda activate viewer_env

python view.py

call conda deactivate

pause
