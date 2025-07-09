import sys
import os
import numpy as np

# force params
m = 50
k = 10
ol = 0.25
d = 1

#############################
# main parameters to modify #
#############################

fld = sys.argv[1]
sim = int(sys.argv[2])
# mcParams = np.loadtxt(fld+'/params.csv', delimiter=',')
# params = mcParams[sim, :]
pST = int(sys.argv[3]) #phenotype state transition 
dp = float(sys.argv[4])
kp = float(sys.argv[5])

print("starting simulations")

tcellMigBias = 0.2 #0.076
cd4Diff = 0.3
macMigBias = 0.1
mdscMigBias = 0.1
nkMigBias = 0.125
macM1 = 0.4 # This needs to be fine tuned, but the data from ERT shows more than double M1 than M2
macM2 = 0.2 # This needs to be fine tuned, but the data from ERT shows more than double M1 than M2
cd8RecRate = 0.1# 0.00192 # these are the percentages of T cells from non-tx tumors provided by ERT
cd4RecRate = 0.1# 0.00192 # these are the percentages of T cells from non-tx tumors provided by ERT
macRecRate = 0.1 # 0.00125 # these are the percentages of mature myeloid cells from non-tx tumors provided by ERT
nkRecRate = 0.1 # 0.00044 # these are the percentages of NK cells from non-tx tumors provided by ERT
mdscRecRate = 0.1 # 0.00469 # these are the percentages of MDSCs from non-tx tumors provided by ERT
recDelay = 0.00
necroticGrowth = 0 #params[12]
necrosisLimit = 940.3
nkKillProb = 0.7 # This has been arbitrarily chosen.

anti_pd1_Dose = 1.5
anti_ctla4_dose = 0.6667 # based on conversion from mg to nM of the experimental schedules from R-T group. See documentation

#############################
# ------------------------- #
#############################

cellParams = np.zeros((14, 6))

# cancer params
cellParams[0, 0] = m  # mu
cellParams[1, 0] = k  # kc
cellParams[2, 0] = d  # damping
cellParams[3, 0] = ol  # overlap
cellParams[4, 0] = 20.0  # diameter (um)
cellParams[5, 0] = 1/24.0  # div probability (hours)
cellParams[6, 0] = 1/(24.0*3.5) # death probability (hours)
cellParams[7, 0] = 40.0  # influence distance

# cd4 params
cellParams[0, 1] = m  # mu
cellParams[1, 1] = k  # kc
cellParams[2, 1] = d  # damping
cellParams[3, 1] = ol  # overlap
cellParams[4, 1] = 10.0  # diameter (um)
cellParams[5, 1] = 1/(24.0*3.0) # lifespan (days) # Gong 2017 (for CD8 cells)
cellParams[6, 1] = 200.0  # migration speed base um/hr
cellParams[7, 1] = cd4Diff  # differentiation to Treg
cellParams[8, 1] = 40.0  # influence distance
cellParams[9, 1] = tcellMigBias  # migration bias base

# cd8 params
cellParams[0, 2] = m  # mu
cellParams[1, 2] = k  # kc
cellParams[2, 2] = d  # damping
cellParams[3, 2] = ol  # overlap
cellParams[4, 2] = 10.0  # diameter (um)
cellParams[5, 2] = 1/(24.0*14.0) # 1/lifespan (days) https://pmc.ncbi.nlm.nih.gov/articles/PMC4489929/
cellParams[6, 2] = 200  # migration speed um/hr | https://onlinelibrary.wiley.com/doi/epdf/10.1038/icb.2012.75 -> scaled based on model scale
cellParams[7, 2] = 0.8 # killProb
cellParams[8, 2] = 2.0  # infScale -> arbitrarily set
cellParams[9, 2] = 40.0  # influence distance
cellParams[10, 2] = tcellMigBias  # migration bias base
cellParams[11, 2] = 0.25 # only for testing #0.053  # proliferation prob
cellParams[12, 2] = 2.0 # arbitrary death scale
cellParams[13, 2] = 2.0 # arbitrary migrate scale

# macrophage params
cellParams[0, 3] = m  # mu
cellParams[1, 3] = k  # kc
cellParams[2, 3] = d  # damping
cellParams[3, 3] = ol  # overlap
cellParams[4, 3] = 20.0  # diameter (um)
cellParams[5, 3] = 1/(24.0*3.0) # (set based on cd8) lifespan (days) -> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4040600/ | Petty and Yang, Tumor-associated macrophages: implications in cancer immunotherapy, 2017
cellParams[6, 3] = 200.0  # migration speed um/hr | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702045/ -> scaled based on model scale
cellParams[7, 3] = macM1  # kM1
cellParams[8, 3] = macM2  # kM2
cellParams[9, 3] = 40.0  # influence distance
cellParams[10, 3] = macMigBias  # migration bias


# nk params
cellParams[0, 4] = m  # mu
cellParams[1, 4] = k  # kc
cellParams[2, 4] = d  # damping
cellParams[3, 4] = ol  # overlap
cellParams[4, 4] = 10.0  # diameter (um) NK cells are actually 6-7 microns: https://pmc.ncbi.nlm.nih.gov/articles/PMC4566154/
cellParams[5, 4] = 1/(24.0*14.0) # lifespan (days) -> 10.1111/j.1365-2567.2007.02573.x
cellParams[6, 4] = 200.0  # migration speed um/hr | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702045/ -> scaled based on model scale
cellParams[7, 4] = nkKillProb  # killProb
cellParams[8, 4] = 2  # infScale -> arbitrarily set
cellParams[9, 4] = 40.0  # influence distance
cellParams[10, 4] = nkMigBias  # migration bias
cellParams[11, 4] = 0.0031 # proliferation probability: Lutz_2011_JImmunol
cellParams[12, 4] = 2.0 # arbitrary
cellParams[13, 4] = 2.0 # arbitrary


# mdsc params
cellParams[0, 5] = m  # mu
cellParams[1, 5] = k  # kc
cellParams[2, 5] = d  # damping
cellParams[3, 5] = ol  # overlap
cellParams[4, 5] = 10.0  # diameter (um) https://pubmed.ncbi.nlm.nih.gov/32038621/
cellParams[5, 5] = 1/(24.0*2.0)  # lifespan (days) https://pmc.ncbi.nlm.nih.gov/articles/PMC5765878/
cellParams[6, 5] = 200.0 # migration speed um/hr (assumed to be the same as other cells)
cellParams[7, 5] = 40  # influenceDistance (assumed to be same as cancer cells, can change)
cellParams[8, 5] = mdscMigBias #


recParams = np.zeros((7, 1))
recParams[0] = cd8RecRate # cd8RecRate
recParams[1] = macRecRate # mRecRate
recParams[2] = cd4RecRate # cd4RecRate
recParams[3] = nkRecRate # nkRecRate
recParams[4] = mdscRecRate # mdscRecRate

recParams[5] = 100.0  # recDist (recruit a uniform distribution recDist away from the tumor edge
recParams[6] = recDelay # recruitment delay (days)

envParams = np.zeros((7, 1))
envParams[0] = 15.0  # initTumorSize x | circle radius
envParams[1] = 25.0#5.0 #5.0 # simulation duration (days)
envParams[2] = necroticGrowth # necrotic growth
envParams[3] = 0#0.5 # necrotic region outward force
envParams[4] = necrosisLimit # necrosis limit (accounts for diffusion limit of oxygen, but is adjustable based on the scale of the simulation)
envParams[5] = anti_pd1_Dose
envParams[6] = anti_ctla4_dose # experimental data

# os.system('mkdir -p ' + sys.argv[1] + '/params')
saveFld = sys.argv[1]+'/set_'+sys.argv[2]+'/params'
os.system('mkdir -p '+saveFld)
np.savetxt(saveFld+'/cellParams.csv', cellParams, delimiter=',')
np.savetxt(saveFld+'/recParams.csv', recParams, delimiter=',')
np.savetxt(saveFld+'/envParams.csv', envParams, delimiter=',')

# np.savetxt(sys.argv[1] + '/params/cellParams.csv', cellParams, delimiter=',')
# np.savetxt(sys.argv[1] + '/params/recParams.csv', recParams, delimiter=',')
# np.savetxt(sys.argv[1] + '/params/envParams.csv', envParams, delimiter=',')


print('------Parameters Grid Done-----\n')
