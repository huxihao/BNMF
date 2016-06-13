## this script is used to publish the minimum sources for reproducing
import time, os, shutil

base = '.'
pub = 'public'
archive = 'bnmf%s'%time.strftime('%Y%m%d')

os.system('rm -rf %s/*'%pub)
os.mkdir('%s/data'%pub)
os.mkdir('%s/work'%pub)
os.mkdir('%s/src'%pub)

shutil.copy(base+'/read_me.txt', pub+'/read_me.txt')

files = [
    'tools.py',
    'Science2009SuppCodeIII.py', ## library of theoretical curves
    'ContactMap.py', ## library of map decomposation
    'decompose_maps.py', ## decompose different hi-c data
    'human_maps.py', ## focus on human hi-c data
    ]

for filename in os.listdir(base+'/src'):
#    if filename in files:
    shutil.copy(base+'/src/'+filename, pub+'/src/'+filename)

for filename in os.listdir(base+'/src'):
    if filename.startswith('table'):
        shutil.copy(base+'result/'+filename, pub+'result/'+filename)

## package
shutil.make_archive(archive, 'gztar', base, pub)

print 'Source codes are saved to', archive

if raw_input('Do you want to test the main modules? y/n ') == 'n':
    exit()
os.system('python src/tools.py WorkPath=%s/work'%pub)
print 'End of Testing!'

