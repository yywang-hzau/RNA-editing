import os
import pandas as pd
import re


def dataframe():
    df1=pd.read_csv('ya-transcriptcount.txt',sep=' ',header=None)

    column1 = df1.columns.values.tolist()
    df2=pd.DataFrame()

    for i1 in column1[1:]:
        df2[i1]=df1[i1].apply(lambda x: 0 if x<5 else 1)
    df2['Col_sum'] = df2.apply(lambda x: x.sum(), axis=1)
    df2['id']=df1[0]
    df2=df2[df2['Col_sum']>1]

    df3=pd.DataFrame()
    column2 = df2.columns.values.tolist()

    for i2 in column2[0:4]:
        df3[i2]=df2[i2].apply(lambda x: '' if x==0 else 'at_mapping_mpec'+str(i2)+'.bam')
    df3['id']=df2['id']

    df3.to_csv('1.txt',sep='\t',index=None,header=None)

def trans_to_sashimi():
    with open ('1.txt','r') as f1:
        list1=f1.readlines()

    with open('2.txt','w') as f2:
        for i in list1:
            i=i.strip()
            seq=i.split('MYZPE13164_O_EIv2.1')[0].strip()
            seq1=','.join(re.split(r'\s+',seq))
            seq2='MYZPE13164_O_EIv2.1'+i.split('MYZPE13164_O_EIv2.1')[1]
            f2.write(seq1+'\t'+seq2+'\n')


def trans_gff():
    with open('2.txt','r') as f1:
        list1=f1.readlines()
    with open('stringtie_merged.gff3','r') as f2:
        list2=f2.readlines()
    with open('3.txt','w') as f3:
        for i in list1:
            i=i.strip()
            for i2 in list2:
                i2=i2.strip()
                if 'transcript' in i2:
                    if i.split('\t')[1] in i2:
                        pos=i2.split('\t')[0]+':'+i2.split('\t')[6]+':'+i2.split('\t')[3]+':'+i2.split('\t')[4]
                        f3.write(i+'\t'+pos+'\n')



def main():
    os.chdir('/Users/yuanyuan/Downloads/cyz/sashimi/')
    dataframe()
    trans_to_sashimi()
    trans_gff()
main()