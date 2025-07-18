import os
import pandas as pd
import portion as P

def merge_bed(dir):
    df1 = pd.read_csv('merge.txt',sep='\t')

    for file in os.listdir(dir):
        if os.path.isfile(file) and file.endswith('bed') :
            df2=pd.read_csv(file,sep='\t',header=None)
            df2.columns=['chr','pos-1','pos',file.split('099')[0],file+'pvalue','strand']
            df2.drop(['pos-1'],axis=1, inplace=True)
            df2.drop(['strand'],axis=1, inplace=True)
            df2=df2[df2[file+'pvalue']>0.99]
            df2.drop([file + 'pvalue'], axis=1, inplace=True)
            df1=pd.merge(df1, df2, on=['chr','pos'], how='outer')
    df1.to_csv('merge_bed.txt',sep='\t',index=None)

def gtf_to_region():
    exon = {}
    cds = {}
    list_transcript = []

    with open('mpec.gtf', 'r') as f1:
        #有#号开头的行跳过
        for line1 in f1:
            if line1.isspace() or line1.startswith('#'):
                continue
            #创建字典，把确定的cds exon放进去
            else:
                list_cds = []
                list_exon = []
                id = line1.split('\t')[0] + '\t' + \
                     line1.split('transcript_id')[1].split(' gene_id')[0].replace(';','').replace('"', '').strip() + \
                     '\t' + line1.split('\t')[6]

                #列表里放的是transcript的区间
                if 'transcript' + '\t' in line1:
                    list_transcript.append(id + '|' + line1.split('\t')[3] + ',' + line1.split('\t')[4])
                elif 'exon' in line1:
                    list_exon.append(line1.split('\t')[3] + ',' + line1.split('\t')[4])
                elif 'CDS' in line1:
                    list_cds.append(line1.split('\t')[3] + ',' + line1.split('\t')[4])

            if id not in exon:
                exon[id] = list_exon
            else:
                exon[id] += list_exon
            if id not in cds:
                cds[id] = list_cds
            else:
                cds[id] += list_cds
    #以染色体+pos为循环的依据，chrom.txt放的是染色体编号
    with open ('chrom.txt','r') as chr_f:
        for chr_id in chr_f.readlines():
            #每条染色体，使区间清零
            chr_id=chr_id.strip()
            exon_extra=P.empty()
            intron=P.empty()
            three=P.empty()
            five=P.empty()
            cds2=P.empty()
            #每条转录本分好exon cds 3`utr 5`utr intron并放入区间里
            for t in list_transcript:
                if chr_id == t.split('\t')[0]:
                    id = str(t.split('|')[0])
            # transcript区间，inter0是转录本
                    t1 = int(t.split('|')[1].split(',')[0])
                    t2 = int(t.split('|')[1].split(',')[1])
                    inter0 = P.closed(t1, t2)
            # exon的第一个区间,inter2是外显子
                    n1 = int(exon[id][0].split(',')[0])

                    n2 = int(exon[id][0].split(',')[1])
                    inter2 = P.closed(n1, n2)
                    for seq_exon in exon[id][1:]:
                        seq_1 = int(seq_exon.split(',')[0])
                        seq_2 = int(seq_exon.split(',')[1])
                        inter2 = inter2 | P.closed(seq_1, seq_2)
            # cds如果为空，直接外显子
                    if cds[id] == []:
                        print(id + '\t' + 'no cds')
                        exon_extra = exon_extra | inter2
                    #cds如果不为空再计算utr
                    else:
                        i1 = int(cds[id][0].split(',')[0])
                        i2 = int(cds[id][0].split(',')[1])
                        inter1 = P.closed(i1, i2)
                        for seq in cds[id][1:]:
                            seq1 = int(seq.split(',')[0])
                            seq2 = int(seq.split(',')[1])
                            inter1 = inter1 | P.closed(seq1, seq2)
                        cds2=cds2 | inter1

                #utr是外显子-CDS
                        utr = inter2 - inter1
                        if '+' in id:
                            five_utr = utr[0]
                            three_utr = utr[-1]
                        elif '-' in id:
                            five_utr = utr[-1]
                            three_utr = utr[0]
                        three=three | three_utr
                        five= five | five_utr
                #内含子是转录本-外显子
                        intron = intron | (inter0 - inter2)
                #除utr外其他的
                        exon_extra = exon_extra | (utr - five_utr - three_utr)

            #开始判断位点是否在区间内
            with open('region.txt', 'a+') as f5:
                #这里没有列名，文件生成后需要自己加上'chr'+'\t'+'pos'+'\t'+'genomic_region',再做剩下的
                with open('merge_bed.txt', 'r') as f2:
                #写入最终文件，只保留了染色体，位置，还有区间信息，可用于后续merge
                    list2 = f2.readlines()
                    #这里跳过第一行，因为是pos表头
                    for i in list2[1:]:
                        i=i.strip()
                        #染色体必须是==,因为如果是in的话，111/11/1都会出现
                        if chr_id == i.split('\t')[0]:
                            pos = i.split('\t')[1]
                            if int(pos) in exon_extra:
                                print(chr_id+'\t'+pos + '\t' + 'exon',file=f5)
                            elif int(pos) in intron:
                                print(chr_id+'\t'+pos+'\t'+'intron',file=f5)
                            elif int(pos) in three:
                                print(chr_id+'\t'+pos+'\t'+'3`UTR',file=f5)
                            elif int(pos) in five:
                                print(chr_id+'\t'+pos+'\t'+'5`UTR',file=f5)
                            elif int(pos) in cds2:
                                print(chr_id+'\t'+pos+'\t'+'CDS',file=f5)

def merge_pos_region():
    df1=pd.read_csv('region.txt',sep='\t',low_memory=False)
    df2=pd.read_csv('merge_bed.txt',sep='\t',low_memory=False)
    df3=pd.merge(df1,df2,how='right',on=['chr','pos'])
    df3.to_csv('merge_bed_region.txt',sep='\t')


def count_region():
    #带有region的进行分组求和
    df1=pd.read_csv('merge_bed_region.txt', sep='\t',low_memory=False)
    df1 = df1.drop_duplicates()
    df2 = pd.DataFrame(df1.groupby(df1['genomic_region']).count().T)
    df2=df2.reset_index()
    #求sum
    df4 = pd.DataFrame(df1.count())
    df4=df4.reset_index()
    df4.columns=['index','sum']
    #merge
    df5=pd.merge(df2,df4,how='outer',on='index')
    df5['ratio_3']=df5['3`UTR']/df5['sum']
    df5['ratio_5'] = df5['5`UTR'] / df5['sum']
    df5['ratio_cds'] = df5['CDS'] / df5['sum']
    df5['ratio_exon'] = df5['exon'] / df5['sum']
    df5['ratio_intron'] = df5['intron'] / df5['sum']

    df5['sum_ratio']=df5['ratio_3']+df5['ratio_5']+df5['ratio_cds']+df5['ratio_exon']+df5['ratio_intron']
    df5['ratio_intergenic'] = 1 - df5['sum_ratio']
    df5.to_csv('merge_bed_region_count.txt',sep='\t')

def merge_bed_099_nonsyn(dir):
    f1 = open('count_syn.txt', 'a+')
    f1.write('nonsyn_count'+'\t'+'sum'+'\t'+'nonsyn_ratio'+'\t'+'syn_ratio'+'\t'+'sample'+'\n')
    for file in os.listdir(dir):
        if os.path.isfile(file) and file.endswith('099.bed'):
            name = file.split('mp')[0]
            df1 = pd.read_csv(file, sep='\t', header=None)
            df1.columns = ['chr', 'pos-1', 'pos', name.split('mp')[0], name + 'pvalue', 'strand']
            df1.drop(['pos-1','strand',name+'pvalue'], axis=1, inplace=True)

            df2 = pd.read_csv(name + 'ann_nonsyn.vcf', sep=' ', header=None)
            df2.columns = ['chr', 'pos']
            df3 = pd.merge(df1, df2, on=['chr', 'pos'], how='inner')
            df3.to_csv(name + 'non_syn_099.txt', sep='\t', index=None)
            nonsyn_ratio=len(df3)/len(df1)
            syn_ratio=1-nonsyn_ratio
            print(str(len(df3)) + '\t'+str(len(df1)) + '\t' + str(nonsyn_ratio)+'\t'+str(syn_ratio)+'\t'+name, file=f1)
    f1.close()


def merge_upsetr(dir):
    df1 = pd.read_csv('merge.txt', sep='\t')
    for file in os.listdir(dir):
        if os.path.isfile(file) and file.endswith('non_syn_099.txt'):
            df2 = pd.read_csv(file, sep='\t')
            df1 = pd.merge(df1, df2, on=['chr', 'pos'], how='outer')
            df1 = df1.fillna(0)
            column1 = df1.columns.values.tolist()
            for i1 in column1[2:]:
                df1[i1] = df1[i1].apply(lambda x: 0 if x == 0 else 1)
    #df1.to_csv('upsetr.txt', sep='\t', index=None)

    df3 = pd.read_csv('merge.txt', sep='\t')
    for file in os.listdir(dir):
        if os.path.isfile(file) and file.endswith('non_syn_099.txt'):
            df4 = pd.read_csv(file, sep='\t')
            df3 = pd.merge(df3, df4, on=['chr', 'pos'], how='outer')
            df3 = df3.fillna(0)
            column1 = df3.columns.values.tolist()
            for i1 in column1[2:]:
                df3[i1] = df3[i1].apply(lambda x: str(x.split('|')[-1]) if x != 0 else 0)
    #df3.to_csv('merge_editinglevel.txt', sep='\t', index=None)
    for i2 in column1[2:]:
        df3=df3.drop(df3[df3[i2]==0].index)

    df3.to_csv('no_zero.txt',sep='\t',index=None)

def change_title():
    df1 = pd.read_csv('rnaseqlist.txt', sep='\t')
    df2 = pd.read_csv('count_removedna.txt', sep='\t')
    df2 = df2.set_index(keys=['group'])
    df3 = pd.DataFrame(df2.T)
    df3['id'] = df3.index
    df4 = pd.merge(df1, df3, how='inner', on='id')
    #df4 = df4.groupby(['group']).sum()
    #column1 = df4.columns.values.tolist()
    #for i1 in column1[1:]:
        #df4[i1] = df4[i1].apply(lambda x: 0 if x == 0 else 1)
    #df6 = df4.T
    df4.to_csv('count_rnasnp_removedna_list.txt', sep='\t',index=None)


def reflect_to_genes():
    with open('merge_editinglevel.txt','r') as f1:
        list1=f1.readlines()
    with open('genes_n.gff','r') as f4:
        list4=f4.readlines()
    with open('nonsyn_gene_editing_level.txt','a+') as f11:
        for i in list1[1:]:
            i=i.strip()
            chr=i.split('\t')[0]
            pos=i.split('\t')[1]
            for i4 in list4:
                i4=i4.strip()
                left=i4.split(' ')[1]
                right=i4.split(' ')[2]

                if chr in i4:
                    if int(pos) > int(left) and int (pos) < int(right):
                        id=i4.split('geneID=')[1]
                        f11.write(i+'\t'+id+'\n')


def reorder_col():
    df1=pd.read_csv('/Users/yuanyuan/Downloads/cyz/RNA-editing/tpm_result/tpm_gt_77.txt',sep='\t',index_col=0)
    df2=pd.read_csv('/Users/yuanyuan/Downloads/cyz/RNA-editing/tpm_result/tpm_ag_52.txt',sep='\t',index_col=0)
    df3=pd.read_csv('/Users/yuanyuan/Downloads/cyz/RNA-editing/tpm_result/tpm_ct_151.txt',sep='\t',index_col=0)
    with open('/Users/yuanyuan/Downloads/cyz/RNA-editing/tpm_result/list-tpm.txt','r') as f1:
        list1=f1.readlines()
    list2=[]
    for i in list1:
        i=i.strip()
        list2.append(i)
    df1=df1[list2]
    df2=df2[list2]
    df3=df3[list2]
    df1.to_csv('/Users/yuanyuan/Downloads/cyz/RNA-editing/tpm_result/tpm_gt_77_final.txt',sep='\t')
    df2.to_csv('/Users/yuanyuan/Downloads/cyz/RNA-editing/tpm_result/tpm_ag_52_final.txt',sep='\t')
    df3.to_csv('/Users/yuanyuan/Downloads/cyz/RNA-editing/tpm_result/tpm_ct_151_final.txt',sep='\t')

def main():
    #修改主路径
    dir = '/Users/yuanyuan/Downloads/cyz/RNA-editing/mpmpn/'
    os.chdir(dir)
    #merge_bed(dir)
    #gtf_to_region()
    #merge_pos_region()
    #count_region()
    #merge_bed_099_nonsyn(dir)
    #merge_upsetr(dir)
    #change_title这一步适用于分组很多的RNAseq
    #change_title()
    reflect_to_genes()
    #reorder_col()

main()