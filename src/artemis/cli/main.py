import os
import subprocess
import sys

def main():
    this_dir = os.path.dirname(__file__)
    package_root = os.path.abspath(os.path.join(this_dir, '..'))  # go up from cli/
    exe_path = os.path.join(package_root, 'bin', 'artemis_executable')
    subprocess.run([exe_path] + sys.argv[1:])