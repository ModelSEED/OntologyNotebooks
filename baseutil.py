import platform
import sys
import sys
import json
from json import dump
import os
import re
from os.path import exists
from pathlib import Path
import logging
import shutil

print("python version " + platform.python_version())

from configparser import ConfigParser

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

config = ConfigParser()
if not exists(str(Path.home()) + '/.kbase/config'):    
    if exists("/scratch/shared/code/sharedconfig.cfg"):
        shutil.copyfile("/scratch/shared/code/sharedconfig.cfg",str(Path.home()) + '/.kbase/config')
    else:
        print("You much create a config file in ~/.kbase/config before running this notebook. See instructions: https://docs.google.com/document/d/1fQ6iS_uaaZKbjWtw1MgzqilklttIibNO9XIIJWgxWKo/edit")
        sys.exit(1)
config.read(str(Path.home()) + '/.kbase/config')
paths = config.get("DevEnv","syspaths").split(";")
codebase = config.get("DevEnv","codebase",fallback="")
for i,filepath in enumerate(paths):
    if filepath[0:1] != "/":
        paths[i] = codebase+"/"+filepath
sys.path = paths + sys.path

from chenry_utility_module.kbdevutils import KBDevUtils

class BaseUtil:
    def __init__(self):
        self.kbdevutil = None
        self.msrecon = None
        self.annoapi = None
    
    def get_kbdevutil(self,name):
        if not self.kbdevutil:
            self.kbdevutil = KBDevUtils(name) 
        return self.kbdevutil
    
    def get_msrecon(self):
        if not self.msrecon:
            self.msrecon = self.get_kbdevutil().msseedrecon()
        return self.msrecon
    
    def get_annoapi(self):
        if not self.annoapi:
            self.annoapi = self.get_kbdevutil().anno_client(native_python_api=True)
        return self.annoapi