'''
COM5345 Cooperative Communications and Networking
Final Project
111064502, 曾雋卿
'''
import math
import random
import numpy as np
import matplotlib.pyplot as plt

random.seed(1024)
N = 20      # number of UE
R = 500     # cell radius
lmb = 1 / 3 # interarrival rate (unit: 1 / timesolt)
T = 10000   # total timeslots for the simulation
K = 3       # K = 1 + (# of relays per sector)
g = 0       # g = relay power (dBm) - UE power (dBm)
# T should be multiple of K
# height: BS 25m; Relay 10m; UE 1.5m
# carrier frequency fc = 700 MHZ
# breakpoint distance = 4 * fc / c * (h_UE - 1)(h_BS - 1)
# dBP_SD = 24 * 0.5 * 2.33494867 * 4    # from UE to BS 
dBP_SD = 112.077536 
# dBP_SR = 9 * 0.5 * 2.33494867 * 4     # from UE to Relay
dBP_SR = 42.029076
# dBP_RD = 2017.39565                   # from Relay to BS
# For packet error probability computation
r = 0.75 * 2                # equivalent code rate
beta = 0.5 / math.pi * math.sqrt(128 / (4 ** r - 1))
phi = 2 ** r - 1 - 0.5 / beta
psi = 2 ** r - 1 + 0.5 / beta
Pn = -174 + 40 + 10 * math.log10(18)    # noise power in dBm
Pu = 23                                 # UE power in dBm
SNR_0 = Pu - Pn

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
        self.packet = []
    def printS(self):
        print("{0}-th UE in the {1}-th sector:".format(
            int(self.order), self.sector))
        print("x =", "{:.2f}".format(self.x),
              "\ty =", "{:.2f}".format(self.y))
    def PacketGenerate(self):
    # as the order-th UE, begin at (K * order)-th timeslot
        # time starting to wait for packet
        t = K * (self.order) - K * Sector[self.sector]
        while t < T:
            t += random.expovariate(lmb)
            self.packet.append(t)
    # discard packets that arrive before the first transmission
        while self.packet[1] <= K * (self.order):
            self.packet.pop(0)
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
                PLoS = (18 / DR_2d[i] + (1 - 18 / DR_2d[i])
                         * math.exp(DR_2d[i] / (-63)))
                if random.random() < PLoS:
                    self.LoS_R[i] = True
                else:
                    self.LoS_R[i] = False
        # LOS probability to BS
        if DBS_2d <= 18:
            self.LoS_BS = True
        else:
            PLoS = (18 / DBS_2d
            + (1 - 18 / DBS_2d) * math.exp(DBS_2d / (-63)))
            if random.random() < PLoS:
                self.LoS_BS = True
            else:
                self.LoS_BS = False
        # Pathloss (in dB) to relay
        # fc is measuered in GHz => 0.7 GHz
        self.LossR = np.zeros([K - 1])
        for i in range(K - 1):
            # Pathloss given Line of Sight
            if DR_2d[i] < dBP_SR:
                self.LossR[i] = (28 
                + 22 * math.log10(DR_3d[i]) - 3.0980392)
            else:
                self.LossR[i] = (28 + 40 * math.log10(DR_3d[i])
                - 3.0980392 - 29.380583)
            # average pathloss must be nonnegative
            if self.LossR[i] < 0:
                self.LossR[i] = 0
            # w/out LOS, compute NLOS version and compare
            if self.LoS_R[i] == False:
                Loss_NLOS = (148.893292
                + 40.32 * (math.log10(DR_3d[i]) - 3))
                if Loss_NLOS > self.LossR[i]:
                    self.LossR[i] = Loss_NLOS
        # Pathloss (in dB) to BS
        self.LossBS = 0.0
        # Pathloss given Line of Sight
        if DBS_2d < dBP_SD:
            self.LossBS = (28 
            + 22 * math.log10(DBS_3d) - 3.0980392)
        else:
            self.LossBS = (28 + 40 * math.log10(DBS_3d)
            - 3.0980392 - 37.059504)
        # pathloss must be nonnegative
        if self.LossBS < 0:
            self.LossBS = 0
        # w/out LOS, compute NLOS version and compare
        if self.LoS_BS == False:
            Loss_NLOS =  (127.705816
            + 39.086386 * (math.log10(DBS_3d) - 3))
            if Loss_NLOS > self.LossBS:
                self.LossBS = Loss_NLOS
    def InitAoI(self):
        # init with -1 (indicating undeclared)
        self.AOI = np.full([T], -1) 
        # suppose last packet succeed, which
        # at (K * (self.order) - K * Sector[self.sector])
        # depart from UE, so BS decode and update AOI at 
        # K * (self.order) - K * Sector[self.sector] + K
        self.AOI[0] = K * (Sector[self.sector] - self.order - 1)

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
PLr = np.zeros([6, K - 1])  # Pathloss from relay to BS
theta = math.pi / 6
match K:
    case 2:
        for i in range(6):
            # relay 0: at center
            Xr[i][0] = 0.6 * R * math.cos(theta)
            Yr[i][0] = 0.6 * R * math.sin(theta)
            # increment theta to go to the next sector
            theta += math.pi / 3
    case 3:
        for i in range(6):
            # relay 0: 10 degree clockwise from center
            Xr[i][0] = 0.6 * R * math.cos(theta - math.pi / 18)
            Yr[i][0] = 0.6 * R * math.sin(theta - math.pi / 18)
            # relay 1: 10 degree counterclockwise from center
            Xr[i][1] = 0.6 * R * math.cos(theta + math.pi / 18)
            Yr[i][1] = 0.6 * R * math.sin(theta + math.pi / 18)
            theta += math.pi / 3
    case 4: 
        for i in range(6):
            # relay 0: 10 degree clockwise from center
            Xr[i][0] = 0.6 * R * math.cos(theta - math.pi / 18)
            Yr[i][0] = 0.6 * R * math.sin(theta - math.pi / 18)
            # relay 1: 10 degree counterclockwise from center
            Xr[i][1] = 0.6 * R * math.cos(theta + math.pi / 18)
            Yr[i][1] = 0.6 * R * math.sin(theta + math.pi / 18)
            # relay 2: at center
            Xr[i][2] = 0.6 * R * math.cos(theta)
            Yr[i][2] = 0.6 * R * math.sin(theta)
            theta += math.pi / 3
# compute pathloss from relay to BS
# Shadowing is neglected
# P_LOS = 18 / 300 + (1 - 18 / 300) * math.exp(300 / (-63))
P_LOS = 0.06803635
for i in range(6):
    for j in range(K - 1):
        # 0.6 R < breakpoint distance
        # (28 + 22 * log10(sqrt(300 ** 2
        # + 15 ** 2)) - 3.0980392)
        if random.random() < P_LOS:
            PLr[i][j] = 79.4105566
        else:                       # non-LOS
            PLr[i][j] = 102.1895676
    # PLr[i][j] = (161.04 + 0.4 * log10(20)
    # - (24.37 - 2.368) * log10(25)
    # + (43.42 - 3.1 * log10(25)) * (log10(sqrt(90225)) - 3)
    # - 3.0980392 - (3.2 * log10(17.625) * log10(17.625) - 4.97)
    # - 0.6 * 8.5)
S = []                  # list of all the N's UE
# how many UE each sector contains
Sector = np.zeros([6], dtype = int)
for i in range(N):
    S.append(source(X[i], Y[i]))
for i in range(N):      # put UE to corresponding sector
    S[i].order = Sector[S[i].sector]
    Sector[S[i].sector] += 1
for i in range(N):
    S[i].PacketGenerate()   # pre-generate packets
    S[i].PLOS()             # compute pathloss
    S[i].InitAoI()          # initialize AoI record
# a table of indices of UEs in each sector 
SectorUE = np.full([6, np.max(Sector)], -1, dtype = int)
# build this table by checking each UE
J = np.zeros([6], dtype = int)  # indices for sectors
for i in range(N):
    j = S[i].sector     # sector to which this UE belongs 
    SectorUE[j][J[j]] = i
    J[j] += 1           # index increment
aveAoI = 0                  # average AoI
# relays that decode successful
helper = np.zeros([K - 1], dtype = 'b') 
# simulate transmission 
for i in range(6):              # for each sector
    t = 0                       # starting from time 0
    j = 0                       # UE to transmit
    while (t < T - 1):  # take 1 timeslot to reach BS
        u = SectorUE[i][j]      # the UE to transmit
        # 0. determine which packet to send
        # discard outdated packets
        while (S[u].packet[0] < t):
            S[u].packet.pop(0)
        # then send packet[0]
        t += 1
        # 1. broadcast
        # compute SNR received by BS
        alpha = (math.sqrt(random.gauss(0, 1) ** 2
                         + random.gauss(0, 1) ** 2)
                    / math.sqrt(2))
        # subtracted by pathloss
        SNR_dB = SNR_0 - S[u].LossBS
        SNR = alpha * (10 ** (SNR_dB / 10))
        # determine whether BS decode successfully
        ErrProb = (1 - beta * SNR
                       * (math.exp(-1 * phi / SNR)
                        - math.exp(-1 * psi / SNR)))
        if random.random() > ErrProb:       # success
            S[u].AOI[t] = (t - S[u].packet[0])
        # if BS decode successfully, relays can rest
        else:
        # 2. relays transmit
        # compute PER from UE to relays
            for k in range(K - 1):
                # generate Rayleigh fading coefficient
                alpha = (math.sqrt(random.gauss(0, 1) ** 2
                                 + random.gauss(0, 1) ** 2)
                        / math.sqrt(2))
                # subtracted by pathloss
                SNR_dB = SNR_0 - S[u].LossR[j]
                SNR_r = alpha * (10 ** (SNR_dB / 10))
                # PER at relays
                ErrProb = (1 - beta * SNR_r
                           * (math.exp(-1 * phi / SNR)
                            - math.exp(-1 * psi / SNR)))
                # relay succeed to decode => help forward
                if random.random() > ErrProb:
                    helper[j] = True
                else:
                    helper[j] = False
        # relays that decode successfully help
            for j in range(K - 1):
                if helper[j] == True:
                    alpha = (math.sqrt(random.gauss(0, 1) ** 2
                                     + random.gauss(0, 1) ** 2)
                            / math.sqrt(2))
                    # subtracted by pathloss
                    SNR_dB = SNR_0 + g - PLr[S[u].sector][j]
                    # BS perform MRC
                    SNR += alpha * (10 ** (SNR_dB / 10))
                    t += 1
        # BS try to decode again
            ErrProb = (1 - beta * SNR
                       * (math.exp(-1 * phi / SNR)
                        - math.exp(-1 * psi / SNR)))
            if random.random() > ErrProb:       # success
                S[u].AOI[t] = (t - S[u].packet[0])
        S[u].packet.pop(0)      # this packet is sent
        # 3. increment of UE index
        j += 1                  # go to next UE
        if (j >= Sector[i]):
            j = 0