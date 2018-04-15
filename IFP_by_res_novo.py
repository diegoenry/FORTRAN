
# coding: utf-8

# In[8]:


from matplotlib import pyplot as plt
import pandas as pd


# In[9]:


list="Mol|AA30|AA31|AG32|AW33|AY36|AP38|AT39|AQ40|AK41|AW42|AG43|AP44|AQ45|AG46|AD47|AR48|AE49|AH50|AP51|AF61|AG63|AA64|AN78|AT79|AF112|DA151|DA152|DG153|DW154|DS156|DY157|DP159|DT160|DQ161|DK162|DW163|DG164|DP165|DQ166|DG167|DD168|DR169|DE170|DH171|DP172|DF182|DG184|DA185|DN199|DT200|DF233".split("|")


# In[10]:


# Le o arquivo
df = pd.read_csv('ifp.csv',header=None)

# Coloca o nome dos residuos
df.columns=list

# Remove Residues with ZERO interactions
for col in list:
    if df[col].max()==0:
        del df[col]
#del data['DF233']

# Transforma a tabela e plota.
df2=df.transpose()
df2.columns=['Polar','Hydrophobic']
df2[1:].plot(kind='bar')
plt.show()

