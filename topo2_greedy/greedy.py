import time
import copy
import sys
import numpy as np
import itertools
import scipy.io as sio

# global variables
numSimus = 1
numTrafficType = 4
stepLen = 1  # --> 50
# demandMaxVal = 400
maxStaN = 10  # max number of stations

log = open("log.txt", "wa")
topo = 30

maxfileName = "../demandMaxPool.mat"
print "reading demandmax from", maxfileName

mat = sio.loadmat(maxfileName)
demandMaxActual_all = mat['demandMaxPool']

minfileName = "../demandMinPool.mat"
print "reading demandMinActual from", minfileName

mat = sio.loadmat(minfileName)
demandMinActual_all = mat['demandMinPool']

demandMaxActual_all = np.reshape(demandMaxActual_all[-2000:],(2000,4,3))
print "dmdmax shape", np.shape(demandMaxActual_all)
demandMinActual_all = np.reshape(demandMinActual_all[-2000:],(2000,4,3))
print "dmdmin shape", np.shape(demandMinActual_all)
#demandMaxActual_all = demandMinVal + np.asarray([list(itertools.product(range(0, 800, stepLen), repeat=4))])
#demandMaxActual_all = demandMaxActual_all[0].reshape((len(demandMaxActual_all[0]), numTrafficType, 1))

iterable = xrange(len(demandMaxActual_all))
numSimulations = len(demandMaxActual_all)

print "numSimulations", len(demandMaxActual_all)
log.write("numSimulations %d" % len(demandMaxActual_all))

assert np.all(demandMinActual_all <= demandMaxActual_all), "demandMinActual should not be greater than any demandMax"

linkCapaFile = "topo" + str(topo) +"Capa.txt"
topoInfoFile = "topo" + str(topo) +".txt"

# utility functions
# initial values to start with, to be extended to array/matrix

lambda1 = 0.00133;
beta1 = 0;
lambda2 = 0.08;
beta2 = 350;
lambda3 = 0.03651;
beta3 = 0.5;
lambda4 = 0.00229;
beta4 = 1;


def utility_1(x, lambda1, beta1): return x * lambda1 + beta1


def utility_2(x, lambda2, beta2): return 1 / (1 + np.exp(-(x - beta2) * lambda2))


def utility_3(x, lambda3, beta3): return x ** beta3 * lambda3


def utility_4(x, lambda4, beta4): return np.log(x * lambda4 + beta4)


def utility_All(x, lambda_all, beta_all):
    return np.array([utility_1(x[0], lambda_all[0], beta_all[0]), utility_2(x[1], lambda_all[1], beta_all[1]), \
                     utility_3(x[2], lambda_all[2], beta_all[2]), utility_4(x[3], lambda_all[3], beta_all[3])])


# read links capacities from .txt file
with open(linkCapaFile) as f1:
    capacity = f1.readlines()
    capacity = [map(int, x.strip().split("\t")) for x in capacity]

capacity = np.asarray(capacity) / 10e5  # mbps

print "capacity Value: \n", capacity
log.write("capacity Value \n %s \n" % str(capacity))

# read flow paths
with open(topoInfoFile) as f2:
    infoTopo = f2.readlines()
    infoTopo = [(x.strip().split("\t")) for x in infoTopo]

numPath = int(infoTopo[0][1])
numStation = int(infoTopo[0][0])

path = infoTopo[numStation + 2: len(infoTopo)]
path = [map(int, x[:-1]) for x in path[::2]]
path = np.asarray(path)

print "numPath ", numPath, "numSta ", numStation
print "path Value: \n", path
log.write("numPath %d numSta %d \n path Value: \n %s \n" % (numPath, numStation, str(path)))

assert numStation == len(capacity), " Dimensions of the two input files should agree "
assert numPath == len(path), " Check topo*.txt file, number of path should agree "

# Get links from path
path_cp = copy.deepcopy(path)
link = np.asarray(path_cp[0][0:2])

for i in path_cp:
    for j in i[1:-1]:
        i.insert(i.index(j), j)
    link = np.append(link, i)

link = np.reshape(link, (len(link) / 2, 2))
link = np.unique(link, axis=0)
print "link Value: \n", link
log.write('link Value: \n %s \n' % str(link))

# Get if link i is on path j
pathInNode = np.zeros((numStation, len(path)))
pathOutNode = np.zeros((numStation, len(path)))

# and the corresponding 'inNodeCapacity'
pathInNodeCapacity = np.full((numStation, len(path)), np.inf)
pathOutNodeCapacity = np.full((numStation, len(path)), np.inf)

for i in range(len(link)):
    for j in range(len(path)):
        for k in range(len(path[j]) - 1):
            if (link[i][0] == path[j][k]) and (link[i][1] == path[j][k + 1]):
                pathOutNode[link[i][0]][j] = 1
                pathInNode[link[i][1]][j] = 1
                pathOutNodeCapacity[link[i][0]][j] = capacity[link[i][0]][link[i][1]]
                pathInNodeCapacity[link[i][1]][j] = capacity[link[i][0]][link[i][1]]

# Path into Node
print "pathInNode \n", pathInNode
print "pathOutNode \n", pathOutNode

print "pathInNodeCapa \n", pathInNodeCapacity
print "pathOutNodeCapa \n", pathOutNodeCapacity

log.write('pathInNode \n %s pathOutNode \n %s pathInNodeCapa \n %s pathOutNodeCapa \n %s \n'
          % (str(pathInNode), str(pathOutNode), str(pathInNodeCapacity), str(pathOutNodeCapacity)))


# Generate output topology information matrix
capacity2w = np.zeros((maxStaN, maxStaN))
capacity2w[:len(capacity), :len(capacity)] = capacity
print "Output \"capacity.npy\" \n", capacity2w
log.write("Output \"capacity.npy\" \n %s \n" % str(capacity2w))

path2w = np.ones((maxStaN, maxStaN)) * (-1.0)
for i in range(len(path)): path2w[i][:len(path[i])] = path[i]
print "Output \"path.npy\" \n", path2w
log.write("Output \"path.npy\" \n %s \n" % str(path2w))

np.save("capacity.npy", capacity2w)
np.save("path.npy", path2w)


# Flows are splitted into those for tenant 1 and for tenant 2
# tenant1_admi is a matrix with one-to-one correspondence to flow rate of different traffic types and different paths
# denoting whether the flow belongs to the tenant1 or not
# tenant1_admi = np.array([[1,0,1],[1,0,1],[1,0,0],[0,1,0]])
#tenant1_admi = np.random.randint(2, size=(numTrafficType, numPath))
tenant1_admi = np.array([[1,0,1],[1,0,1],[1,0,0],[0,1,0]])
np.save("tenant1_admittedFlows", tenant1_admi)
tenant2_admi = 1 - tenant1_admi

print "tenant 1 admitted flows\n", tenant1_admi
print "tenant 2 admitted flows\n", tenant2_admi
log.write("tenant 1 admitted flows\n %s \n tenant 2 admitted flows\n %s \n" % (str(tenant1_admi), str(tenant2_admi)))

# parameters for the 2 tenants are same for now (different across traffic types),
# can be adjusted in the future
lambda_t1 = tenant1_admi * np.array([[lambda1], [lambda2], [lambda3], [lambda4]])
lambda_t2 = tenant2_admi * np.array([[lambda1], [lambda2], [lambda3], [lambda4]])
# lambda_t2 = tenant2_admi * np.array([[lambda1 +1], [lambda2+1], [lambda3+1], [lambda4+1]])
beta_t1 = tenant1_admi * np.array([[beta1], [beta2], [beta3], [beta4]])
beta_t2 = tenant2_admi * np.array([[beta1], [beta2], [beta3], [beta4]])

# Sum lambda_1 and 2 for the ease of utility computation
lambda_all = lambda_t1 + lambda_t2
beta_all = beta_t1 + beta_t2
print "lambda for tenant 1:\n", lambda_t1
print "lambda_all\n", lambda_all

log.write("lambda for tenant 1:\n %s \n lambda_all\n %s \n" % (str(lambda_t1), str(lambda_all)))


utilityMax = [0]
rateLog = np.zeros((1, numTrafficType, numPath))
demandMaxLog = np.zeros((1, numTrafficType, numPath))
demandMinLog = np.zeros((1, numTrafficType, numPath))
timeUsed = [0]
# less relavent but interesting to see
inTimeMax = np.zeros((1, numStation))
outTimeMax = np.zeros((1, numStation))

print len(demandMaxActual_all)

for simIter in range(len(demandMaxActual_all)):
    # preparation for serious looping
    rateMaxUtil = np.zeros((numTrafficType, numPath))
    inNodeTimeSumMax = np.zeros((1, numStation))
    outNodeTimeSumMax = np.zeros((1, numStation))

    # upper demand in this loop
    print "/*******start simulation ", simIter, "*********/"
    demandMaxActual = demandMaxActual_all[simIter]
    demandMinActual = demandMinActual_all[simIter]
    
    print "Traffic type specific demand \n", demandMaxActual
    
    print "start searching...."

    timeStart = time.clock()

    maxFound = False
    rate = np.ones((numTrafficType, numPath)) * demandMinActual
    print "rate0: ",rate
    
    utilityTotal = 0
    utility_current = 0
    # flows whose rates not yet constrained
    increAllowed = np.ones((numTrafficType, numPath))

    while maxFound == False and np.any(np.nonzero(increAllowed)):
        # potential increment in utility if increase a single rate
        poIncreaseU = utility_All(rate + stepLen, lambda_all, beta_all) - utility_All(rate, lambda_all, beta_all)
        poIncreaseU = poIncreaseU * increAllowed

        [typeId, pathId] = np.argwhere(poIncreaseU == np.amax(poIncreaseU))[0]

        poRate = copy.deepcopy(rate)

        poRate[typeId][pathId] = poRate[typeId][pathId] + stepLen

        # check time constraint first
        pathSumRate = np.sum(poRate, axis=(0))  # each element shall not be greater than the capacity (linkCapacity[])
        inNodeTimeSum = pathInNode * pathSumRate / pathInNodeCapacity
        inNodeTimeSum = np.sum(inNodeTimeSum, axis=1)
        outNodeTimeSum = pathOutNode * pathSumRate / pathOutNodeCapacity
        outNodeTimeSum = np.sum(outNodeTimeSum, axis=1)

        # if time is fine
        if (np.all(np.less_equal(inNodeTimeSum, np.ones(np.shape(inNodeTimeSum))))) and (
                np.all(np.less_equal(outNodeTimeSum, np.ones(np.shape(outNodeTimeSum))))):
            # confirm increase
            rate[typeId][pathId] = min(poRate[typeId][pathId], demandMaxActual[typeId][pathId])

            if rate[typeId][pathId] == demandMaxActual[typeId][pathId]:
                increAllowed[typeId][pathId] = 0

        else:
            increAllowed[typeId][pathId] = 0
            continue

        assert np.all(np.less_equal(rate, demandMaxActual)), " some demand exceed"

        util1 = utility_1(rate[0], lambda_all[0], beta_all[0])
        # print util1
        util2 = utility_2(rate[1], lambda_all[1], beta_all[1])
        # print util2
        util3 = utility_3(rate[2], lambda_all[2], beta_all[2])
        # print util3
        util4 = utility_4(rate[3], lambda_all[3], beta_all[3])
        # print util4
        utilityTotal = np.sum([util1 + util2 + util3 + util4])
        # compute utility and check if increase from previous
#        utilityTotal = np.sum(utility_All(rate, lambda_all, beta_all))

        if utilityTotal > utility_current:
            utility_current = utilityTotal
            rateMaxUtil = copy.deepcopy(rate)
            inNodeTimeSumMax = copy.deepcopy(inNodeTimeSum)
            outNodeTimeSumMax = copy.deepcopy(outNodeTimeSum)

            #        elif np.any(np.nonzero(increAllowed)) and utilityTotal <= utility_current:
            #            maxFound = True
            # exit(-1), "utility not increasing"

    # counter.value += 1
    # if counter.value > 0 and counter.value % int(len(demandMaxActual_all) / 10.0) == 0:
    #     print ('Simulation checkpoint %d, %d%% completed' % (
    #     counter.value, counter.value / int(len(demandMaxActual_all) / 100.0)))
    #     log.write('Simulation checkpoint %d, %s%% completed \n' % (
    #     counter.value, str(counter.value / int(len(demandMaxActual_all) / 100.0))))

    rateLog[-1] = rateMaxUtil
    rateLog = np.append(rateLog, np.zeros((1, numTrafficType, numPath)), axis=0)
    demandMaxLog[-1] = demandMaxActual
    demandMaxLog = np.append(demandMaxLog, np.zeros((1, numTrafficType, numPath)), axis=0)
    demandMinLog[-1] = demandMinActual
    demandMinLog = np.append(demandMinLog, np.zeros((1, numTrafficType, numPath)), axis=0)

    utilityMax[-1] = utility_current
    utilityMax = np.append(utilityMax,[0],axis=0)
    timeUsed[-1] = time.clock() - timeStart
    timeUsed = np.append(timeUsed,[0])

    inTimeMax[-1] = inNodeTimeSumMax
    inTimeMax = np.append(inTimeMax,np.zeros((1, numStation)), axis=0)
    outTimeMax[-1] = outNodeTimeSumMax
    outTimeMax = np.append(outTimeMax, np.zeros((1, numStation)), axis=0)

#    rateLog_update = rateLog.get()
#    rateLog_update[-1] = rateMaxUtil
#    rateLog_update = np.append(rateLog_update, np.zeros((1, numTrafficType, numPath)), axis=0)
#    rateLog.put(rateLog_update)
#
#    demandMaxLog_update = demandMaxLog.get()
#    demandMaxLog_update[-1] = demandMaxActual
#    demandMaxLog_update = np.append(demandMaxLog_update, np.zeros((1, numTrafficType, numPath)), axis=0)
#    demandMaxLog.put(demandMaxLog_update)
#
#    demandMinLog_update = demandMinLog.get()
#    demandMinLog_update[-1] = demandMinActual
#    demandMinLog_update = np.append(demandMinLog_update, np.zeros((1, numTrafficType, numPath)), axis=0)
#    demandMinLog.put(demandMinLog_update)
#
#    utilityMax_update = utilityMax.get()
#    utilityMax_update[-1] = utility_current
#    utilityMax_update = np.append(utilityMax_update, [0], axis=0)
#    utilityMax.put(utilityMax_update)
#
#    timeUsed_update = timeUsed.get()
#    timeUsed_update[-1] = time.clock() - timeStart
#    timeUsed_update = np.append(timeUsed_update, [0])
#    timeUsed.put(timeUsed_update)
#
#    inTimeMax_update = inTimeMax.get()
#    inTimeMax_update[-1] = inNodeTimeSumMax
#    inTimeMax_update = np.append(inTimeMax_update, np.zeros((1, numStation)), axis=0)
#    inTimeMax.put(inTimeMax_update)
#    outTimeMax_update = outTimeMax.get()
#    outTimeMax_update[-1] = outNodeTimeSumMax
#    outTimeMax_update = np.append(outTimeMax_update, np.zeros((1, numStation)), axis=0)
#    outTimeMax.put(outTimeMax_update)

    np.save("data/utilityTotal", utilityMax[:-1])
    np.save("data/demandMaxLog", demandMaxLog[:-1])
    np.save("data/demandMinLog", demandMinLog[:-1])
    np.save("data/rateLog", rateLog[:-1])
    np.save("data/time", timeUsed[:-1])
    np.save("data/inTimeMax", inTimeMax[:-1])
    np.save("data/outTimeMax", outTimeMax[:-1])

    print "time used: ", timeUsed[-2], "sec"
    print "optimal utility ", utility_current, " at rate \n", rateMaxUtil
