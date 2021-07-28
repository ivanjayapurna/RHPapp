#libraries
import csv
import pandas as pd
import numpy as np
import math
from scipy.optimize import fsolve

#exp_name (name it), [ mon1, mon2, mon3, mon4, ] (monomer names)
# [ mon1_mr, mon2_mr, mon3_mr, mon4_mr, ] (monomer proportions out of 100)
#int_std (name of int_std), int_std_mr (ratio of int_std), 
#init (name of initiator), init_conc (concentration of initiator), 
#raft (name of raft), macro_raft_mw (molecular proportion of raft), i_r_ratio (unknown, probably 1), 
#solv (name of solvent), conversion (conv), target_mw, 
#rxn_scale (prolly should be 1), rxn_vol (L), rxn_temp (C)


def main(exp_name, mon1, mon2, mon3, mon4, mon1_mr, mon2_mr, mon3_mr, mon4_mr, int_std, int_std_mr, init, init_conc, raft, macro_raft_mw, i_r_ratio, solv, conversion, target_mw, rxn_scale, rxn_vol, rxn_temp):
	###############
	#  FUNCTIONS  #
	###############

	# helper function used to replace sum() from excel
	def calc_total(arr, r):
		total = 0
		for i in range(r):
			total += arr[i]
		return total

	# main body function to calculate inital values for output reactor scheme
	def initial_values():
		mols[0] = (mon1_mr/100)*(rxn_scale/1000) #mol
		mols[1] = (mon2_mr/100)*(rxn_scale/1000)
		mols[2] = (mon3_mr/100)*(rxn_scale/1000)
		mols[3] = (mon4_mr/100)*(rxn_scale/1000)
		mols[4] = (int_std_mr/100)*(rxn_scale/1000)
		mols[5] = init_mr*calc_total(mols,5)/(100 + int_std_mr)
		mols[6] = raft_mr*calc_total(mols,5)/(100 + int_std_mr)
		for i in range(7):
			if reagents[i] != '':
				try:
					if (reagents[i] == 'MACRORAFT'):
						MWs[i] = macro_raft_mw #g/mol
						concs[i] = mols[i]/(rxn_vol/1000) #mol/L
						weights[i] = mols[i]*MWs[i]*1000 #mg
						vols[i] = weights[i]/(liq_ps[i]*1000) #mL
					else:
						MWs[i] = df.loc[reagents[i],'M (g/mol)']
						if df.loc[reagents[i],'p (g/mL)'] < 100:
							liq_ps[i] = df.loc[reagents[i],'p (g/mL)']
						concs[i] = mols[i]/(rxn_vol/1000) #mol/L
						weights[i] = mols[i]*MWs[i]*1000 #mg
						if i == 5:
							vols[i] = weights[i]/(init_conc) #mL = mg/(mg/mL)
						else:
							vols[i] = weights[i]/(liq_ps[i]*1000) #mL = mg/(mg/mL)
				except KeyError:
					print('ERROR: Make sure all input reagent codes are valid.')


	# Simultaneous equation system to be solved via fsolve
	# this essentially replicates excel goal seek
	def f_sys(x):
		# define variables
		_raft_mr = x[0]
		_init_mr = x[1]
		_curr_mw = x[2]
		_mols5 = x[3]
		_mols6 = x[4]
		_concs0 = x[5]
		_concs1 = x[6]
		_concs2 = x[7]
		_concs3 = x[8]
		_concs4 = x[9]
		_concs5 = x[10]
		_concs6 = x[11]
		_weights0 = x[12]
		_weights1 = x[13]
		_weights2 = x[14]
		_weights3 = x[15]
		_weights4 = x[16]
		_weights5 = x[17]
		_weights6 = x[18]
		_vols0 = x[19]
		_vols1 = x[20]
		_vols2 = x[21]
		_vols3 = x[22]
		_vols4 = x[23]
		_vols5 = x[24]
		_vols6 = x[25]

		# define functions
		return [_init_mr - _raft_mr/i_r_ratio,
				_mols5 - _init_mr*(mols[0] + mols[1] + mols[2] + mols[3] + mols[4])/(100 + int_std_mr), #mol
				_mols6 - _raft_mr*(mols[0] + mols[1] + mols[2] + mols[3] + mols[4])/(100 + int_std_mr), #mol
				_concs0 - mols[0]*1000/rxn_vol, #mol/L
				_concs1 - mols[1]*1000/rxn_vol, #mol/L
				_concs2 - mols[2]*1000/rxn_vol, #mol/L
				_concs3 - mols[3]*1000/rxn_vol, #mol/L
				_concs4 - mols[4]*1000/rxn_vol, #mol/L
				_concs5 - _mols5*1000/rxn_vol, #mol/L
				_concs6 - _mols6*1000/rxn_vol, #mol/L
				_weights0 - mols[0]*MWs[0]*1000, #mg
				_weights1 - mols[1]*MWs[1]*1000, #mg
				_weights2 - mols[2]*MWs[2]*1000, #mg
				_weights3 - mols[3]*MWs[3]*1000, #mg
				_weights4 - mols[4]*MWs[4]*1000, #mg
				_weights5 - _mols5*MWs[5]*1000, #mg
				_weights6 - _mols6*MWs[6]*1000, #mg
				_vols0 - _weights0/(liq_ps[0]*1000), #mL = mg/((g/mL)*(mg/g))
				_vols1 - _weights1/(liq_ps[1]*1000), #mL
				_vols2 - _weights2/(liq_ps[2]*1000), #mL
				_vols3 - _weights3/(liq_ps[3]*1000), #mL
				_vols4 - _weights4/(liq_ps[4]*1000), #mL
				_vols5 - _weights5/(init_conc), #mL = mg/(mg/mL)
				_vols6 - _weights6/(liq_ps[6]*1000), #mL
				_curr_mw - ((_weights0+_weights1+_weights2+_weights3+_weights4)*conversion)/(_mols6*1000) - MWs[6],
				_curr_mw - target_mw ] #mg

	# reassign values to output arrays based on simultaneous eq solver results
	def reassign(x):
		raft_mr = x[0]
		init_mr = x[1]
		mols[5] = x[3]
		mols[6] = x[4]
		concs[0] = x[5]
		concs[1] = x[6]
		concs[2] = x[7]
		concs[3] = x[8]
		concs[4] = x[9]
		concs[5] = x[10]
		concs[6] = x[11]
		weights[0] = x[12]
		weights[1] = x[13]
		weights[2] = x[14]
		weights[3] = x[15]
		weights[4] = x[16]
		weights[5] = x[17]
		weights[6] = x[18]
		vols[0] = x[19]
		vols[1] = x[20]
		vols[2] = x[21]
		vols[3] = x[22]
		vols[4] = x[23]
		vols[5] = x[24]
		vols[6] = x[25]

	############
	#  SCRIPT  #
	############

	# open up reagent master list as a new pandas data frame
	df = pd.read_csv('data/reagent_master_list.txt', delimiter='\t')
	df.index = df['Code']
	
	# first check input monomer molecular ratios sum to 100. if not print error statement
	if (mon1_mr + mon2_mr + mon3_mr + mon4_mr) != 100:
		print ('ERROR: Make sure monomer molecular ratios sum to 100.')
	else:
		# calculate initial / tentative values for all reagents (not solvent or total)
		raft_mr = 0.5
		init_mr = raft_mr/i_r_ratio
		reagents = [mon1, mon2, mon3, mon4, int_std, init, raft, solv, 'Total']
		mols = [0, 0, 0, 0, 0, 0, 0, 0, 0]
		MWs = [0, 0, 0, 0, 0, 0, 0, 0, 0]
		liq_ps = [float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), float('inf')]
		weights = [0, 0, 0, 0, 0, 0, 0, 0, 0]
		vols = [0, 0, 0, 0, 0, 0, 0, 0, 0]
		concs = [0, 0, 0, 0, 0, 0, 0, 0, 0]
		initial_values()
		# run fsolve to calculate values properly via solving simultaneous, based on target_mw
		curr_mw = (calc_total(weights,5)*conversion)/(mols[6]*1000) + MWs[6]
		print('initial guess mw: ' + str(curr_mw)) # FOR TESTING, CAN DELETE LATER
		x0 = [raft_mr, init_mr, curr_mw, mols[5], mols[6], concs[0], concs[1], concs[2], concs[3],
				concs[4], concs[5], concs[6], weights[0], weights[1], weights[2], weights[3], weights[4],
				weights[5], weights[6], vols[0], vols[1], vols[2], vols[3], vols[4], vols[5], vols[6]]
		x = fsolve(f_sys, x0, factor=0.1)
		reassign(x)
		
		curr_mw = (calc_total(weights,5)*conversion)/(mols[6]*1000) + MWs[6]

		# calculate the solvent and total column properties 
		vols[7] = (rxn_vol - calc_total(vols, 7)) #mL
		weights[7] = (vols[7]*df.loc[reagents[7],'p (g/mL)']*1000) #mg
		mols[7] = (weights[7]/(df.loc[reagents[7],'M (g/mol)']*1000)) #mol
		concs[7] = (mols[7]/(rxn_vol/1000)) #mol
		vols[8] = (calc_total(vols, 8)) #mL
		weights[8] = (calc_total(weights, 8)) #mg
		concs[8] = (calc_total(concs, 8)) #M (i.e. mol/L)
		mols[8] = (calc_total(mols, 8)) #mol

		# organize output dataframe
		d = {'vol (mL)':vols, 'mass (mg)':weights, 'conc (M)':concs, 'mols (mol)':mols}
		df2 = pd.DataFrame(data=d)
		pd.set_option('display.expand_frame_repr', False)
		r_names = [df.loc[x, 'Name'] if x != 'Total' else 'Total' for x in reagents]
		df2.index = r_names

		dp = 100*(target_mw - MWs[6])/(mon1_mr*MWs[0] + mon2_mr*MWs[1] + mon3_mr*MWs[2] + mon4_mr*MWs[3] + int_std_mr*MWs[4])
		
		# temporary output of reaction time
		# constants
		init_df = pd.read_csv('data/initiators.txt', delimiter='\t')
		init_df.index = init_df['Code']
		kd = np.exp(init_df.loc[init,'Freq Factor (ln A)'])*np.exp(-init_df.loc[init,'Ea (KJ/mol K)']*1000/(8.3144621*(rxn_temp+273)))
		kp = 10e3
		ka = 10e6
		kf = 10e4
		kt = 10e7
		kct = 10e7
		T0 = concs[6]
		I0 = concs[5]*2.0 # MULTIPLY BY 2 OR NAH??
		M0 = 0
		for i in range(4):
			M0 += abs(concs[i])
		print('kd: ', kd)
		print('T0: ', T0)
		print('I0: ', I0)
		print('M0: ', M0)
		# lumped parameters
		alpha = 1.0 + kt*kf/(2.0*T0*kct*ka)
		beta = 2.0*alpha*alpha*T0/I0 - 1.0
		gamma = kp*np.sqrt(kf*alpha/(kd*kct*ka))
		delta = ((np.sqrt(1.0 + beta) - 1.0)/(np.sqrt(1.0 + beta) + 1.0)) * ((1.0-conversion)**(-1.0/gamma))
		# finally time itself (Wang paper eq rearranged)
		time = (1.0/kd)*np.log((1.0/beta)*(((1.0+delta)/(1.0-delta))**2.0 - 1.0))
		print('time: ' + str(time))

		# return output values to main webapp
		return exp_name, target_mw, curr_mw, df2.round(4), dp, time
