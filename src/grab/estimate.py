import numpy as np
import pway_funcs as fn2
import scipy.io
import GRAB as grab
import os

max_iter = 1000
K = 20

filenames = os.listdir('../../cache/data/')

for fn in filenames:
  skip = False
  for skipname in ['real', 'erdos-renyi', 'hub']:
    if fn.startswith(skipname):
      skip = True
  if skip:
    continue
  print(fn)
  folder_net = 'net_' + str(K)
  folder_blocks = 'blocks_' + str(K)
  if not os.path.exists('../../cache/grab/output/' + folder_net):
    os.makedirs('../../cache/grab/output/' + folder_net)
  if not os.path.exists('../../cache/grab/output/' + folder_blocks):
    os.makedirs('../../cache/grab/output/' + folder_blocks)

  for lmbda in [x * 0.01 for x in range(10, 50, 2)]:
    print(lmbda)
    if not os.path.isfile('../../cache/grab/output/' + folder_net + '/' + fn + '_' + str(lmbda)):
      S = np.corrcoef(np.loadtxt('../../cache/grab/data/' + fn).T)
      (Theta, blocks) = grab.BCD(S,lmbda=lmbda,K=K, max_iter=max_iter)

      with open('../../cache/grab/output/' + folder_blocks + '/' + fn + '_' + str(lmbda), 'w') as filehandle:  
        for listitem in blocks:
          for innerlistitem in listitem:
            filehandle.write('%s ' % str(innerlistitem))
          filehandle.write('\n')

      net = np.zeros(shape=Theta.shape)
      net[Theta != 0] = 1
      np.savetxt('../../cache/grab/output/' + folder_net + '/' + fn + '_' + str(lmbda), net)
