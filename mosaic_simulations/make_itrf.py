import pandas as pd
import sys

array = sys.argv[1:]

## load arrays
for i in range(len(array)):
	if i == 0:
		df = pd.read_csv('%s.itrf'%array[i],delimiter=" ", header=None,names=['X', 'Y', 'Z', 'dish_diam', 'station', 'mount'],index_col=False)
	else:
		df2 = pd.read_csv('%s.itrf'%array[i],delimiter=" ", header=None,names=['X', 'Y', 'Z', 'dish_diam', 'station', 'mount'],index_col=False)
		df = df.append(df2,ignore_index=True)

df.to_csv('%s_sims.itrf'%"_".join(array),header=False,index=False,sep=' ')

## load vla array

df_vla = pd.read_csv('vlab.itrf',delimiter=" ", header=0,names=['X', 'Y', 'Z', 'dish_diam', 'station', 'mount'],index_col=False)

for i in range(len(df)):
	df['X'][i] = df_vla.iloc[i]['X']
	df['Y'][i] = df_vla.iloc[i]['Y']
	df['Z'][i] = df_vla.iloc[i]['Z']
df.to_csv('%s_vlapos_sims.itrf'%"_".join(array),header=False,index=False,sep=' ')