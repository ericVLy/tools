import os
import subprocess

subprocess.Popen("pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple", shell=True)
subprocess.Popen("pip install -r requirements.txt", shell=True)
