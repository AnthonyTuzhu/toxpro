import os, glob, ntpath
import pandas as pd
from flask import send_from_directory

BASE = os.getenv('TOXPRO_FILE_DIR')


class DB:

    pass