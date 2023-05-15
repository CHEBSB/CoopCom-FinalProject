'''
COM5345 Cooperative Communications and Networking
Final Project
111064502, 曾雋卿
'''
import math
import random
import numpy as np

N = 20      # number of UE
R = 500     # cell radius
lmb = 1 / 3 # interarrival rate (unit: 1 / timesolt)
T = 50      # total timeslots for the simulation
K = 3       # K = 1 + (number of relays per sector)
# T should be multiple of K
# height: BS 25m; Relay 10m; UE 1.5m
# carrier frequency fc = 700 MHZ

class source:
    def __init__(self, X, Y):
        self.order = 0
        self.x = X
        self.y = Y
        # to which of the 6 sectors it belongs
        theta = math.atan2(Y, X)
        if theta >= 0 and theta < math.pi / 3:
            self.sector = 0
        elif theta >= math.pi / 3 and theta < math.pi * 2 / 3:
            self.sector = 1
        elif theta >= math.pi * 2 / 3 and theta <= math.pi:
            self.sector = 2
        elif theta >= -1 * math.pi and theta < -2 * math.pi / 3:
            self.sector = 3
        elif theta >= -2 * math.pi / 3 and theta < math.pi / (-3):
            self.sector = 4
        else:
            self.sector = 5
        # timings when packets are generated
        self.packetTime = []
    def printS(self):
        print("{0}-th UE in the {1}-th sector:".format(
            int(self.order), self.sector))
        print("x =", "{:.3f}".format(self.x),
              "\ty =", "{:.3f}".format(self.y))
    def PacketGenerate(self):
    # as the order-th UE, begin at (K * order)-th timeslot
        # time starting to wait for packet
        t = K * (self.order) - K * Sector[self.sector]
        while t < T:
            t += random.expovariate(lmb)
            self.packetTime.append(t)
    # discard packets that won't be sent
        # related to # of UE in this sector
        step = K * Sector[self.sector]
        # timeslots at which this source can transmit
        TxTime = range(K * int(self.order), T + 1, step)
        i = len(self.packetTime) - 1
        j = len(TxTime) - 1
        # packets arrive after the last transmission
        while self.packetTime[i] > TxTime[j]:
            self.packetTime.pop(i)
            i -= 1
        # timeslot of the last transmission
        while self.packetTime[i] <= TxTime[j - 1]:
            j -= 1
        # not-the-freshest packets
        i -= 1
        while i >= 0 and j >= 1:
            if self.packetTime[i] > TxTime[j - 1]:
                self.packetTime.pop(i)
            else:
                while j >= 1 and self.packetTime[
                    i] <= TxTime[j - 1]:
                    j -= 1
            i -= 1
        # packet arrive before the first transmission
        while self.packetTime[1] <= K * (self.order):
            self.packetTime.pop(0)
    # compute LOS probability and Pathloss
    def PLOS(self):
        # distance to each relay
        DR_2d = np.zeros([K - 1])
        DR_3d = np.zeros([K - 1])
        for i in range(K - 1):
            DR_2d[i] = math.sqrt((self.x - Xr[self.sector][i])**2 
                            + (self.x - Yr[self.sector][i])**2)
            # height diff = 8.5m; 8.5 ** 2 = 72.25
            DR_3d[i] = math.sqrt((DR_2d[i])**2 + 72.25)
        # distance to BS
        DBS_2d = math.sqrt((self.x)**2 + (self.y)**2)
        DBS_3d = math.sqrt((DBS_2d)**2 + 552.25)
        # LOS probability from this source to relay
        self.LoS_R = np.zeros([K - 1], dtype = 'b')
        for i in range(K - 1):
            if DR_2d[i] <= 18:
                self.LoS_R[i] = True
            else:
                PLoS = (18 / DR_2d[i]
                + (1 - 18 / DR_2d[i]) * math.exp(DR_2d[i] / (-63)))
                if random.random() < PLoS:
                    self.LoS_R[i] = True
                else:
                    self.LoS_R[i] = False
        # LOS probability to BS
        if DR_2d[i] <= 18:
            self.LoS_BS = True
        else:
            PLoS = (18 / DBS_2d
            + (1 - 18 / DBS_2d) * math.exp(DBS_2d / (-63)))
            if random.random() < PLoS:
                self.LoS_BS = True
            else:
                self.LoS_BS = False
        # Pathloss to relay
            

# Locations of UE
X = np.array([-278.4520829528112, -267.652186218928,
-166.44373105820478, 32.301111191102336, 265.3104022644901,
265.626537368792, 147.4183645398856, -159.67072898278036,
-129.2403955435645, -108.60532626016919, -280.89843872820586,
-3.821438258458727, 30.817742218392368, -245.65608450129938,
-114.10733406682527, 319.17568984958257, -144.89855881070326,
438.4660415849662, 363.2847046544206, -246.88861639109982])
Y = np.array([36.47148887813148, -376.4822964934558,
354.9143976350996, 5.536338423434927, 280.458995497302,
-254.5281714303589, -331.7786948600575, 250.31381971501332,
379.43594387732367, -265.51883471640815, -335.35750663156557,
-117.07630897682532, -295.32968595822604, 59.48309753052422,
-234.46432177001606, 260.6832269364211, -47.958294309139546,
-96.45572615294537, 199.25795360348002, -390.168611516649])
# Location of Relays
Xr = np.zeros([6, K - 1])
Yr = np.zeros([6, K - 1])
theta = math.pi / 6
for i in range(6):
    # relay 0: 10 degree clockwise from center
    Xr[i][0] = 0.6 * R * math.cos(theta - math.pi / 18)
    Yr[i][0] = 0.6 * R * math.sin(theta - math.pi / 18)
    # relay 1: 10 degree counterclockwise from center
    Xr[i][1] = 0.6 * R * math.cos(theta + math.pi / 18)
    Yr[i][1] = 0.6 * R * math.sin(theta + math.pi / 18)
    # if exist: relay 2: at the center
    # Xr[i][2] = 0.6 * R * math.cos(theta)
    # Yr[i][2] = 0.6 * R * math.sin(theta)
    # increment theta to go to the next sector
    theta += math.pi / 3

S = []                  # list of all the N's UE
# how many UE each sector contains
Sector = np.zeros([6], dtype = int) 
for i in range(N):
    S.append(source(X[i], Y[i]))
for i in range(N):      # put UE to corresponding sector
    S[i].order = Sector[S[i].sector]
    Sector[S[i].sector] += 1
    S[i].PacketGenerate()   # pre-generate packets
S[1].PLOS()