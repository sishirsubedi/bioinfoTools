
import pandas as pd

df_indel = pd.read_csv('varscan.somatic.indel.vcf.txt',sep='\t',skiprows=17)
df_snp = pd.read_csv('varscan.somatic.snp.vcf.txt',sep='\t',skiprows=17)

def getInfo(df):
    result = []
    for indx,row in df.iterrows():

      temp=[]
      temp.append(row['#CHROM'])
      temp.append(row['POS'])
      temp.append(row['REF'])
      temp.append(row['ALT'])


      info=row['INFO'].split(';')
      for i in info:
        s_line=i.split('=')
        if s_line[0]=='SS':
          temp.append(s_line[1])
        if s_line[0]=='SPV':
          temp.append(s_line[1])

      for j in row['NORMAL'].split(':'):
        temp.append(j)

      for j in row['TUMOR'].split(':'):
        temp.append(j)

      result.append(temp)

    return result

df_snp_filter=pd.DataFrame(getInfo(df_snp))
df_snp_filter.columns=['chrom','position','ref','var','status','somatic_pval',\
'n_genotype','n_g_quality','n_depth','n_r_depth','n_v_depth','n_v_freq','n_dp4',\
't_genotype','t_g_quality','t_depth','t_r_depth','t_v_depth','t_v_freq','t_dp4']

df_indel_filter=pd.DataFrame(getInfo(df_indel))
df_indel_filter.columns=['chrom','position','ref','var','status','somatic_pval',\
'n_genotype','n_g_quality','n_depth','n_r_depth','n_v_depth','n_v_freq','n_dp4',\
't_genotype','t_g_quality','t_depth','t_r_depth','t_v_depth','t_v_freq','t_dp4']


df_res = pd.concat([df_snp_filter, df_indel_filter], ignore_index=True)
df_res.to_csv("bcm_snp_indel.txt",sep='\t',index=False)

#### atlas from bcm
df_indel = pd.read_csv('atlas-indel.somatic.target.indel.txt',sep='\t')
df_indel = df_indel[['Chromosome','Start_position','Reference_Allele','Tumor_Seq_Allele2']]
df_indel.columns=['chrom','position','ref','var']

df_snp = pd.read_csv('atlassnp.somatic.target.snp.txt',sep='\t')
df_snp = df_snp[['Chromosome','Start_position','Reference_Allele','Tumor_Seq_Allele2']]
df_snp.columns=['chrom','position','ref','var']

df_res = pd.concat([df_snp, df_indel], ignore_index=True)
df_res.to_csv("bcm_snp_indel_atlas.target.txt",sep='\t',index=False)
