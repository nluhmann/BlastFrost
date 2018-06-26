#virtualenv in ~/my_software/BlastFrost/.venv/

import matplotlib as mpl
mpl.use('Agg')


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys

sns.set(style="white")


input_file = sys.argv[1]

df = pd.read_csv(input_file,header=0,sep='\t', dtype={'subspecies': str})



df["query"] = df["gene"] + " " + df["query"]



# split strain column to remove paths
#df['strain'] = df[0].str.split('/').str[-1]
#df.drop(df.columns[0],axis=1, inplace=True)

# set strain column as index
#df = df.set_index('strain')

# sum number of 1's in each row, then sort rows in descending order
#df['sum'] = df[list(df.columns)].sum(axis=1,numeric_only=True)
#df.sort_values(by=['sum'],ascending=False, inplace=True)

#df.drop(['sum'],axis=1,inplace=True)


dims = (8, 150.27)
fig, ax = plt.subplots(figsize=dims)

#data = df.pivot("strain","gene","presence")
#print(data)
subset = df.pivot("strain","query","presence")
x = df['strain']
#subset = subset.set_index(['strain','query'])
#subset.sort_index(x['strain'])
subset = subset.reindex(x)
subset.reindex(sorted(subset.columns), axis=1)


# Draw the heatmap
g = sns.heatmap(subset, cmap="Reds", vmax=1, vmin=0, square=False, cbar=False,linewidths=0.1, linecolor='black')

# Ajust axis labels
g.tick_params(axis='y',labelsize=3,labelrotation=45)
g.tick_params(axis='x',labelsize=3)
g.set_xlabel('queried genes',fontsize=5)
g.set_ylabel('reference strains',fontsize=5)


plt.savefig("./heatmap.pdf", format="pdf")

#plt.show()
