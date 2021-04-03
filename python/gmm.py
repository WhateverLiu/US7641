
import sys
import numpy as np
from sklearn.mixture import GaussianMixture


if 0:
  dataPath = sys.argv[1]
  Nclust = int(sys.argv[2])
  savePath = sys.argv[3]
  Nrepeat = 30
  X = np.genfromtxt(dataPath, delimiter = ',')
  bics = np.zeros(Nrepeat)
  labels = []
  for i in range(Nrepeat):
    gm = GaussianMixture(n_components = Nclust, random_state = i * 3 + 7).fit(X)
    labels.append(gm.predict(X))
    bics[i] = gm.bic(X)
  tmp = np.argmin(bics) 
  print(bics)
  bic = bics[tmp]
  lab = labels[tmp]
  rst = np.zeros(len(lab) + 1)
  rst[0] = bic
  rst[1:] = lab
  np.savetxt(savePath, rst)




if 1:
  trainDataPath = sys.argv[1]
  testDataPath = sys.argv[2]
  Nclust = int(sys.argv[3])
  trainSavePath = sys.argv[4]
  testSavePath = sys.argv[5]
  Nrepeat = 30
  Xtrain = np.genfromtxt(trainDataPath, delimiter = ',')
  Xtest = np.genfromtxt(testDataPath, delimiter = ',')
  bics = np.zeros(Nrepeat)
  gmms = []
  for i in range(Nrepeat):
    gm = GaussianMixture(n_components = Nclust, random_state = i * 3 + 7).fit(Xtrain)
    bics[i] = gm.bic(Xtrain)
    gmms.append(gm)
  tmp = np.argmin(bics)
  gm = gmms[tmp]
  trainScore = gm.predict_proba(Xtrain)
  testScore = gm.predict_proba(Xtest)
  np.savetxt(trainSavePath, trainScore, delimiter = ",")
  np.savetxt(testSavePath, testScore, delimiter = ",")





