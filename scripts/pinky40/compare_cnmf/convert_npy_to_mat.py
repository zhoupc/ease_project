import numpy as np 
import scipy.io 

rlt = np.load('3pass_rlts_for_pc.npy')
scipy.io.savemat('cnmf_scan1_results_from_ding.mat', {'rlt':rlt})
