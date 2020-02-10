import os

adirs = os.listdir('.')

dirsum = 0
for adir in adirs:
  os.chdir(adir)
  dirsum = dirsum + len(os.listdir('.'))
  os.chdir('../')
print(dirsum)
