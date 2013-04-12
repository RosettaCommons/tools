import os
from os.path import basename
import glob

for f in glob.glob('./bin/*'):
    os.remove(f)

for dirpath, dirnames, filenames in os.walk('./'):
    if len( basename(dirpath) ) == 0 or basename(dirpath)[0] == '.' or basename(dirpath) == 'bin':
        continue
    for f in filenames:
        filename = os.path.join(dirpath, f)
        if filename[-3:] == '.py':
            print filename
            os.symlink('../' + filename, './bin/' + f)
