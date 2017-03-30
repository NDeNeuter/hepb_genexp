#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

path = '/Users/nicolasdeneuter/Bestanden/PhD/Projects/GOA/RNAseq/PIMS HBV anon.xlsx'
data = pd.read_excel(path, sheetname='Gegevens deelnemers')
X = data[['Geslacht', 'Leeftijd start studie']].loc[1:41,:]
X['Geslacht'] = X['Geslacht'].apply(lambda x: 1 if x == 'V' else 0)
X['Leeftijd start studie'] = X['Leeftijd start studie']/X['Leeftijd start studie'].max()
sns.clustermap(X)
sns.plt.show()