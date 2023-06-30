import pandas as pd
import sys

antennae = sys.argv[1:]

## load arrays
df = pd.read_csv('simulations/master.itrf',delimiter=" ", header=None,names=['X', 'Y', 'Z', 'dish_diam', 'station', 'mount'],index_col=False)
for i in range(len(antennae)):
	if i == 0:
		df3 = df.loc[(df['station'] == antennae[i])].reset_index(drop=True)
	else:
		df3 = df3.append(df.loc[(df['station'] == antennae[i])].reset_index(drop=True),ignore_index=True)

df3.to_csv('sims.itrf',header=False,index=False,sep=' ')

## load vla array

df_vla = pd.read_csv('simulations/vlab.itrf',delimiter=" ", header=0,names=['X', 'Y', 'Z', 'dish_diam', 'station', 'mount'],index_col=False)

for i in range(len(df3)):
	df3['X'][i] = df_vla.iloc[i]['X']
	df3['Y'][i] = df_vla.iloc[i]['Y']
	df3['Z'][i] = df_vla.iloc[i]['Z']
df3.to_csv('vlapos_sims.itrf',header=False,index=False,sep=' ')