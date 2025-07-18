#author: WANG Yuanyuan

import os
import pandas as pd
import portion as P

#如果没有文件夹，建立文件夹
def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def trans_it(dir_revbed):
    os.chdir(dir_revbed)
    for files in os.listdir(dir_revbed):
        if os.path.isfile(files) and files.endswith('_rev099.bed'):
            print(files)
            with open(files,'r') as f1:
                list1=f1.readlines()
            name=files.split('_rev099.bed')[0]+'rev_transtofwd_099.bed'
            with open(name,'w') as f2:
                for i in list1:
                    seq=i.split(' ')
                    left=seq[3].split('>')[0]
                    if left == 'A':
                        left2 = 'T'
                    elif left == 'G':
                        left2='C'
                    elif left == 'C':
                        left2='G'
                    elif left == 'T':
                        left2 = 'A'
                    right=seq[3].split('>')[1]
                    if right == 'A':
                        right2 = 'T'
                    elif right == 'G':
                        right2='C'
                    elif right == 'C':
                        right2='G'
                    elif right == 'T':
                        right2 = 'A'
                    seq2=left2+'>'+right2
                    f2.write(seq[0]+'\t'+seq[1]+'\t'+seq[2]+'\t'+seq2+'\t'+seq[4]+'\t'+seq[5]+'\n')

def merge_fwdrev(dir_fwdbed,dir_revbed,dir_final):
    os.chdir(dir_revbed)
    for files in os.listdir(dir_revbed):
        if os.path.isfile(files) and files.endswith('rev_transtofwd_099.bed'):
            print(files)
            sample_name = files.split('rev_transtofwd_099.bed')[0]
            df1=pd.read_csv(files,sep='\t',header=None)
            print(df1.head())
            #df1.drop([1,4,5], axis=1, inplace=True)
            #fwd文件夹，加上样品名，加上后缀
            fwd = dir_fwdbed + sample_name + '_fwd099.bed'
            df2=pd.read_csv(fwd,sep=' ',header=None)
            print(df2.head())
            #df2.drop([1, 4, 5], axis=1, inplace=True)
            df3=pd.merge(df2,df1,on=[0,1],how='outer')
            df3.columns=['chr','pos','fwd_dp','fwd_type','fwd_editing_level','fwd_p',
                         'rev_dp','rev_type','rev_editing_level','rev_p']
            df3.fillna('/',inplace=True)
            df3['editing_type_fwd_rev']=df3['fwd_type']+'[F]&[R]'+df3['rev_type']
            #df3.drop(['fwd_editing_type','rev_editing_type'], axis=1, inplace=True)
            df3.to_csv(dir_final+sample_name+'_final.bed',index=None,sep='\t')

def bed_editingtype(dir_final):
    os.chdir(dir_final)
    for file in os.listdir(dir_final):
        if os.path.isfile(file) and file.endswith('_final.bed'):
            print(file)
            name=file.split('_final.bed')[0]
            with open(file,'r') as f1:
                list1=f1.readlines()
            with open(name+'_editing_type.bed','w') as f2:
                f2.write('chr'+'\t'+'pos'+'\t'+'editing_level'+'\t'+'editing_type'+'\t'+'strand'+'\n')
                for i in list1[1:]:
                    i=i.strip()
                    chr=i.split('\t')[0]
                    pos=i.split('\t')[1]

                    if '[F]&[R]/' in i:
                        strand='F'
                        type=i.split('\t')[3]
                        level=i.split('\t')[4]
                    elif '/[F]&[R]' in i:
                        strand='R'
                        type = i.split('\t')[7]
                        level = i.split('\t')[8]
                    elif '[F]&[R]' in i and '/' not in i:
                        if float(i.split('\t')[2]) > float(i.split('\t')[6]):
                            strand='Fr'
                            type = i.split('\t')[3]
                            level = i.split('\t')[4]
                        else:
                            strand='Rf'
                            type = i.split('\t')[7]
                            level = i.split('\t')[8]
                    f2.write(chr+'\t'+pos+'\t'+level+'\t'+type+'\t'+strand+'\n')





def merge_bed(dir,dir_final):
    os.chdir(dir_final)
    df1 = pd.read_csv(dir_final+'merge.txt', sep='\t')

    #合并总表格

    for file in os.listdir(dir_final):
        if os.path.isfile(file) and file.endswith('_editing_type.bed'):
            print(file)
            df2 = pd.read_csv(file, sep='\t')
            df2.columns = ['chr', 'pos', file.split('_editing_type.bed')[0]+'level','type', 'strand']

            df1 = pd.merge(df2, df1, on=['chr', 'pos','type','strand'], how='outer')
            df1=df1[(df1['chr']=='scaffold_1')|(df1['chr']=='scaffold_2')|(df1['chr']=='scaffold_3')
                    |(df1['chr']=='scaffold_4')|(df1['chr']=='scaffold_5')|(df1['chr']=='scaffold_6')]

    df1.to_csv(dir+'merge_bed0617.txt', sep='\t', index=None)




def gtf_to_region(dir):
    #切换路径到主目录，下面放:gtf，chrom.txt,merge过的bed文件，这里以上一步生成的merge_bed.txt来做，如果用merge_bed_type.txt也行只是提供位置
    os.chdir(dir)

    #这里的exon放的是没有cds区间的，所以最后命名会改成ncRNA
    exon = {}
    cds = {}
    list_transcript = []

    with open('mpec.gtf', 'r') as f1:
        #有#号开头的行跳过
        for line1 in f1:
            #每一行对应的都是某个基因的信息
            if line1.isspace() or line1.startswith('#'):
                continue
            #创建字典，把确定的cds exon放进去
            else:
                list_cds = []
                list_exon = []

                #把染色体位置、基因名、以及正负链信息提取出来，用于后面判断3'还是5'UTR
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

    #到这一步，cds，exon字典中放了每个基因的区间了，transcript因为是每个基因只有一行，所以放入列表中就行


    #以染色体+pos为循环的依据，chrom.txt放的是染色体编号
    with open ('chrom.txt','r') as chr_f:
        for chr_id in chr_f.readlines():

            #每条染色体，使区间清零，这样做是因为每条染色体的位置区间都重新开始计数
            chr_id=chr_id.strip()
            exon_extra=P.empty()
            ncrna = P.empty()
            intron=P.empty()
            three=P.empty()
            five=P.empty()
            cds2=P.empty()

            #每条转录本分好exon cds 3`utr 5`utr intron并放入区间里
            for t in list_transcript:
                #首先判断转录本是否是这条染色体上的
                if chr_id == t.split('\t')[0]:
                    id = str(t.split('|')[0])

                    #transcript区间，inter0是转录本
                    t1 = int(t.split('|')[1].split(',')[0])
                    t2 = int(t.split('|')[1].split(',')[1])
                    inter0 = P.closed(t1, t2)

                    #inter2是外显子
                    n1 = int(exon[id][0].split(',')[0])
                    n2 = int(exon[id][0].split(',')[1])
                    inter2 = P.closed(n1, n2)

                    for seq_exon in exon[id][1:]:
                        seq_1 = int(seq_exon.split(',')[0])
                        seq_2 = int(seq_exon.split(',')[1])
                        inter2 = inter2 | P.closed(seq_1, seq_2)

                    #cds如果为空，直接定义为外显子(ncRNA)

                    if cds[id] == []:
                        print(id + '\t' + 'no cds')
                        ncrna = ncrna | inter2
                        intron = intron | (inter0 - inter2)

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
            with open('region_726.txt', 'a+') as f5:
                #这里没有列名，文件生成后需要自己加上'chr'+'\t'+'pos'+'\t'+'genomic_region',再做剩下的
                with open('list_726.txt', 'r') as f2:
                #写入最终文件，只保留了染色体，位置，还有区间信息，可用于后续merge
                    list2 = f2.readlines()
                    #这里跳过第一行，因为是pos表头
                    for i in list2[1:]:
                        i=i.strip()
                        #染色体必须是==,因为如果是in的话，111/11/1都会出现
                        if chr_id == i.split('\t')[0]:
                            pos = i.split('\t')[1]
                            if int(pos) in exon_extra:
                                print(chr_id+'\t'+pos + '\t' + 'exon_extra',file=f5)
                            elif int(pos) in intron:
                                print(chr_id+'\t'+pos+'\t'+'intron',file=f5)
                            elif int(pos) in three:
                                print(chr_id+'\t'+pos+'\t'+'3`UTR',file=f5)
                            elif int(pos) in five:
                                print(chr_id+'\t'+pos+'\t'+'5`UTR',file=f5)
                            elif int(pos) in cds2:
                                print(chr_id+'\t'+pos+'\t'+'CDS',file=f5)
                            elif int(pos) in ncrna:
                                print(chr_id + '\t' + pos + '\t' + 'ncrna', file=f5)

def merge_pos_region(dir):
    os.chdir(dir)
    df1 = pd.read_csv('region.txt', sep='\t', low_memory=False)
    df1.columns = ['chr', 'pos', 'region']
    #df1['pos']=df1['pos'].astype('object')
    #同理这里也可以用merge_bed.txt
    df2 = pd.read_csv('merge_bed_type.txt', sep='\t', low_memory=False)

    df3 = pd.merge(df1, df2, how='right', on=['chr','pos'])
    df3.to_csv('merge_bed_region.txt', sep='\t',index=None)


def count_region(dir):
    os.chdir(dir)
    #带有region的进行分组求和
    df1=pd.read_csv('merge_bed_region.txt', sep='\t',low_memory=False,index_col=[0,1])

    df2 = pd.DataFrame(df1.groupby(df1['region']).count().T)
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
    df5['ratio_ncRNA'] = df5['ncrna'] / df5['sum']
    df5['ratio_intron'] = df5['intron'] / df5['sum']

    df5['sum_ratio']=df5['ratio_3']+df5['ratio_5']+df5['ratio_cds']+df5['ratio_ncRNA']+df5['ratio_intron']
    df5['ratio_intergenic'] = 1 - df5['sum_ratio']
    df5.to_csv('merge_bed_region_count.txt',sep='\t',index=None)
    #终文件作图前可能需要简单处理一下，去掉多余的行之类的

def concat_nonsyn(dir_ann,dir):
    os.chdir(dir_ann)
    #ann文件夹下每个txt文件来自于一个vcf文件，因为sailor 是正负链分开的，所有有两个最终的vcf生成,且未过滤conf
    #这里对syn/non_syn的合并不区分样品了，统计在所有样品中出现的nonsyn即可
    list1=[]
    df1=pd.DataFrame()
    for file in os.listdir(dir_ann):
        if os.path.isfile(file) and file.endswith('remove_table.txt'):
            pd_file=pd.read_csv(file,sep='\t',header=None)
            list1.append(pd_file)
            df1=pd.concat(list1)
    df1.columns=['chr','pos','.','ref','var','missense']
    df1.drop(['.','ref','var'], axis=1, inplace=True)
    df1=df1.sort_values(by=['chr','pos'])
    df1=df1.drop_duplicates(subset=['chr','pos'])
    df1.to_csv(dir+'ann.txt',sep='\t',index=None)

def merge_count_nonsyn(dir):
    os.chdir(dir)
    df1=pd.read_csv('ann.txt',sep='\t')
    df1['type']='nonsynonymous'
    df2=pd.read_csv('merge_bed.txt',sep='\t')
    df3=pd.read_csv('merge_bed_type.txt',sep='\t')

    df_ann_bed=pd.merge(df1,df2,how='right',on=['chr','pos'])
    df_ann_bed['type'].fillna('synonymous',inplace=True)
    df_ann_bed.to_csv('merge_bed_ann.txt',sep='\t',index=None)

    df_ann_editing_type=pd.merge(df1,df3,how='right',on=['chr','pos'])
    df_ann_editing_type['type'].fillna('synonymous',inplace=True)
    df_ann_editing_type.to_csv('merge_editing_type_ann.txt',sep='\t',index=None)

    df_ann_bed=df_ann_bed.set_index(['chr','pos','missense'])
    df2 = pd.DataFrame(df_ann_bed.groupby(df_ann_bed['type']).count().T)
    df2.to_csv('count_missense.txt',sep='\t')

def get_editing_level(dir):
    os.chdir(dir)
    df1=pd.read_csv('merge_bed_ann.txt',sep='\t',low_memory=False)

    def f(x):
        if '[F]&[R]/' in x:
            return str(x.split('[F]&[R]/')[0].split('|')[-1])
        elif '/[F]&[R]' in x:
            return str(x.split('|')[-1])
        elif '[F]&[R]' in x:

            left=x.split('[F]&[R]')[0].split('|')[-1]
            right=x.split('|')[-1]

            editing_level=(float(left)+float(right))/2
            return str(editing_level)
        else:
            return 'NA'

    column1 = df1.columns.values.tolist()
    for i1 in column1[4:]:
        df1[i1] = df1[i1].apply(lambda x: f(str(x)))
    df1.to_csv('editing_level.txt',sep='\t',index=None)

def main():
    dir = '/Users/yuanyuan/Downloads/RNA-editing/5hosts/'
    dir_revbed=dir+'rev099bed0617'
    dir_fwdbed=dir+'fwd099bed0617/'
    dir_final=dir+'final_bed0617/'
    dir_ann=dir+'ann/'

    """
    这个脚本的路径需要特别注意一下哦,最好一步一步做，一步一步把#去掉
    主要是整理格式、初步统计结果
    """

    """""
    step1: 准备工作，转换rev文件中的mapping结果为编辑类型,再将一个样品的正负链进行合并，再提取出editing type方便分类、统计
    """""
    #trans_it(dir_revbed)
    #mkdir(dir_final)
    #merge_fwdrev(dir_fwdbed,dir_revbed,dir_final)
    #bed_editingtype(dir_final)

    """""
    step2：合并 bed文件
    """""
    #merge_bed(dir,dir_final)

    """""
    step3: 分类genomic_region，mpec的gtf文件中是没有utr的注释的，所以要人为加上这些区间，
    没有cds的exon注释成ncRNA (这一步速度较慢)
    再将genomic_region和所有的bed文件进行一个合并
    """""
    #gtf_to_region(dir)
    #merge_pos_region(dir)

    """""
    step4:根据region进行统计
    """""
    #count_region(dir)
    """""
    计算好后最后一行是region的作图时去掉
    """""

    """""
    step5:合并所有的nonsynonymous位点
    """""
    #concat_nonsyn(dir_ann,dir)

    """""
    step6:统计突变类型信息
    """""
    #merge_count_nonsyn(dir)

    """""
    step7:RNA editing level提取出来
    """""
    #get_editing_level(dir)

main()
print('finish')