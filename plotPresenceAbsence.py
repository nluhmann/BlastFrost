import matplotlib as mpl
mpl.use('Agg')


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

#plt.use('Agg')

sns.set(style="white")


input_file = sys.argv[1]

df = pd.read_csv(input_file,header=None,sep='\t')

# split strain column to remove paths
df['strains'] = df[0].str.split('/').str[-1]
df.drop(df.columns[0],axis=1, inplace=True)

# set strain column as index
df = df.set_index('strains')

# sum number of 1's in each row, then sort rows in descending order
df['sum'] = df[list(df.columns)].sum(axis=1,numeric_only=True)
df.sort_values(by=['sum'],ascending=False, inplace=True)
df.drop(['sum'],axis=1,inplace=True)

#print(df)

# Draw the heatmap
g = sns.heatmap(df.head(50), cmap="Reds", vmax=1, vmin=0, square=False, linewidths=.5, cbar=False)

# Ajust axis labels
g.tick_params(axis='y',labelsize=3,labelrotation=45)
g.tick_params(axis='x',labelsize=5)
g.set_xlabel('query kmer start position',fontsize=7)
g.set_ylabel('reference strains',fontsize=7)


plt.savefig("./heatmap.pdf", format="pdf")

#plt.show()
