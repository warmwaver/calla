@echo off
cd %cd%
pyinstaller --onefile -i TB_column_64x64.ico TB_column.py
pause