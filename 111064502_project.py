'''
COM5345 Cooperative Communications and Networking
Final Project
111064502, 曾雋卿
'''
import math
import random
import numpy as np

random.seed(1024)
N = 20      # number of UE
R = 500     # cell radius
lmb = 1 / 3 # interarrival rate (unit: 1 / timesolt)
T = 10000   # total timeslots for the simulation
K = 3       # K = 1 + (# of relays per sector)
g = 5       # g = relay power (dBm) - UE power (dBm)
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
Pn = -142 + 40 + 10 * math.log10(18)    # noise power in dBm
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
        # timings for the current packet
        self.packet = -2
        self.nextPacket = 0     # next arrival
        self.AoI_x = []         # timings of AoI update
        self.AoI_y = []         # value for update
    def printS(self):
        print("{0}-th UE in the {1}-th sector:".format(
            int(self.order), self.sector))
        print("x =", "{:.2f}".format(self.x),
              "\ty =", "{:.2f}".format(self.y))
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
    def ComputeAoI(self):
        L = len(self.AoI_x)
        if (L != 0):
            self.AoI = (self.AoI_x[0]) * (self.AoI_x[0])
            for i in range(L - 1):
                self.AoI += ((self.AoI_x[i + 1] - self.AoI_x[i])
                * (self.AoI_y[i] + self.AoI_x[i + 1] - self.AoI_x[i]))
            self.AoI += (T - 1 - self.AoI_x[L - 1]) * (
            self.AoI_y[L - 1] + T - 1 - self.AoI_x[L - 1])
            self.AoI = self.AoI / T / 2 # divide by 2 here
        else:
            self.AoI = (T - 1) * (T - 1) / 2

# Locations of UE
X = np.array([107.22714004313332, -187.28180215443888,
160.92111257559213, 331.7888087894073, -295.2968674814588,
202.24850546797825, -411.7512449735281, -272.0749217659809,
186.29527695369256, 114.7546690567159, -102.09622791055727,
-276.98943749062454, -207.58121318332093, 144.96329871244302,
-360.297097952003, 286.70170211800655, 36.815249114606104,
420.99713196030723, -344.9195875764001, 43.086256651425515])
Y = np.array([-387.950956376743, 385.85275385821797,
-223.4263030453202, -270.5321796955459, -16.79616972940113,
-230.05960111447132, -142.32740629583685, 156.34164822228865,
355.5296500964473, 337.8941549050852, -258.33981178650936,
-333.2136187091145, 343.64085756716497, -54.20884063544719,
198.1077605939539, -6.890779220058789, -3.5768690164863415,
93.16770448698333, -129.716178521889, -249.3647563168273])
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
S = []                  # list of all the N's UE
# how many UE each sector contains
Sector = np.zeros([6], dtype = int)
for i in range(N):
    S.append(source(X[i], Y[i]))
for i in range(N):      # put UE to corresponding sector
    S[i].order = Sector[S[i].sector]
    Sector[S[i].sector] += 1
for i in range(N):
    S[i].PLOS()             # compute pathloss
# a table of indices of UEs in each sector 
SectorUE = np.full([6, np.max(Sector)], -1, dtype = int)
# build this table by checking each UE
J = np.zeros([6], dtype = int)  # indices for sectors
for i in range(N):
    j = S[i].sector     # sector to which this UE belongs 
    SectorUE[j][J[j]] = i
    J[j] += 1           # index increment
# simulate transmission 
for i in range(6):              # for each sector
    t = 0                       # starting from time 0
    j = 0                       # UE to transmit
    while (t < T - K):
        u = SectorUE[i][j]      # the UE to transmit
        # 0. determine whether the packet can be sent
        if S[u].packet > t:     # packet not yet arrive
            t += 1
        else:
            S[u].nextPacket = S[u].packet + random.expovariate(lmb)
            while (S[u].nextPacket <= t):
                S[u].packet = S[u].nextPacket
                S[u].nextPacket += random.expovariate(lmb)
            # then send packet
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
                S[u].AoI_x.append(t)
                S[u].AoI_y.append(t - S[u].packet)
            # if BS decode successfully, relays can rest
            else:
            # 2. relays transmit
            # compute SNR received by relays
                for k in range(K - 1):
                # generate Rayleigh fading coefficient
                    alpha = (math.sqrt(random.gauss(0, 1) ** 2
                                     + random.gauss(0, 1) ** 2)
                            / math.sqrt(2))
                # subtracted by pathloss
                    SNR_dB = SNR_0 - S[u].LossR[k]
                    SNR_r = alpha * (10 ** (SNR_dB / 10))
                # PER at relays
                    ErrProb = (1 - beta * SNR_r
                            * (math.exp(-1 * phi / SNR_r)
                             - math.exp(-1 * psi / SNR_r)))
                # relay decode successfully => help forward
                    if random.random() > ErrProb:
                        alpha = (math.sqrt(random.gauss(0, 1) ** 2
                                         + random.gauss(0, 1) ** 2)
                                / math.sqrt(2))
                        # subtracted by pathloss
                        SNR_dB = SNR_0 + g - PLr[i][k]
                        # BS perform MRC
                        SNR += alpha * (10 ** (SNR_dB / 10))
                        t += 1
        # BS try to decode again
                ErrProb = (1 - beta * SNR
                        * (math.exp(-1 * phi / SNR)
                         - math.exp(-1 * psi / SNR)))
                if random.random() > ErrProb:       # success
                    S[u].AoI_x.append(t)
                    S[u].AoI_y.append(t - S[u].packet)
            S[u].packet = S[u].nextPacket
        # 3. increment of UE index
        j += 1                  # go to next UE
        if (j >= Sector[i]):
            j = 0
aveAoI = 0                  # average AoI
for i in range(N):
    S[i].ComputeAoI()
    aveAoI += S[i].AoI
aveAoI = aveAoI / N
print("average AoI:", "{:.5f}".format(aveAoI))
stdevAoI = 0        # standard deviation
for i in range(N):
    stdevAoI += (S[i].AoI) ** 2
stdevAoI = math.sqrt(stdevAoI / N - aveAoI ** 2)
print("standard dev:", "{:.5f}".format(stdevAoI))