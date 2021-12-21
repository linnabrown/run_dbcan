# resolve path to Hotpep data file,
# both when installed as a package, and in development
# written by Le Nov 14, 2021

import os.path
import pkg_resources

def eCAMI_data_path(*args):
    return pkg_resources.resource_filename('eCAMI', os.path.join(*args))
