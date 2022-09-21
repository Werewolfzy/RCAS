#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Author:Werewolfzy

import sys
import getopt
import pandas as pd
import numpy as np
import os



def usage():
    print(
"""
usage: python [{0}] ... [-g genotype_file | -p phenotype_file | -s significant_snp | -f freq_cutoff | -o out_file]  ...
参数说明:
-g     : 基因型文件bed文件前缀
-p     : 表型文件
-s     : 显著SNP位点信息
-f     : 设定的等位基因频率值，默认是0.3
-o     : 数据保存名字
-h     : 帮助信息
""".format(sys.argv[0]))

def main():

    opts,args = getopt.getopt(sys.argv[1:],"hg:p:s:f:o:")
    genotype_file = ""
    phenotype_file = ""
    significant_snp = ""
    freq_cutoff = 0.3
    out_file = "result.txt"


    for op,value in opts:
        if op == '-g':
            genotype_file = value
        elif op == "-p":
            phenotype_file = value
        elif op == "-s":
            significant_snp = value
        elif op == "-f":
            freq_cutoff = value
        elif op == "-o":
            out_file = value
        else:
            usage()
            sys.exit()
            return


    # print(genotype_file,phenotype_file,significant_snp,freq_cutoff,parent_info,out_file)


    # print(df)

    dfp =  pd.read_csv(phenotype_file,sep='\t')
    dfp.columns = ['ID','type','trait']

    dfp2 = dfp.loc[:,['ID','ID','type','trait']]

    dfp2.columns = ['IID','ID','type','trait']

    dfID = dfp2.loc[:,['IID','ID']]

    dfID.to_csv('test1.txt', sep='\t',index=False, header=None)


    os.system('plink --chr-set 27 --allow-extra-chr --bfile %s --keep test1.txt --recode --make-bed --out test2' %(genotype_file))
    os.system('rm test1.*')
    os.system('plink --chr-set 27 --allow-extra-chr --bfile test2 --extract %s --recode --make-bed --out test3' %(significant_snp))
    os.system('rm test2.*')

    os.system('plink --bfile test3 --recode vcf-iid --out test4')
    os.system('rm test3.*')

    os.system('grep -v "##fileformat*" test4.vcf > test5.vcf')
    os.system('grep -v "##fileDate*" test5.vcf > test4.vcf')
    os.system('grep -v "##source*" test4.vcf > test5.vcf')
    os.system('grep -v "##contig*" test5.vcf > test4.vcf')
    os.system('grep -v "##INFO*" test4.vcf > test5.vcf')
    os.system('grep -v "##FORMAT*" test5.vcf > test4.vcf')
    os.system('cp test4.vcf test5.txt')
    os.system('rm test4.*')


    df5 =  pd.read_csv("test5.txt",sep='\t',header=None)
    df5.drop(columns=[0,1,3,4,5,6,7,8],inplace=True)
    df52 = df5.T
    df52[0] = df52[0].apply(lambda x: x.split('_')[0])
    dfp.to_csv('test6.txt', sep='\t',index=False)
    df6 = pd.read_csv("test6.txt",sep='\t',header=None)
    df7 = pd.merge(df6,df52,how = 'left', on = 0)

    df7.to_csv('test7.txt', sep='\t',index=False, header=None)

    os.system('rm test5.*')
    os.system('rm test6.*')

    df8 =  pd.read_csv("test7.txt",sep='\t')
    df9 = df8.replace('./.','5')


    count1 = (df9['type'] == 1).sum().sum()
    count2 = (df9['type'] == 2).sum().sum()
    resulttxt = df9.iloc[:,0:3]
    resulttxt.rename(columns={'ID':'name'}, inplace = True)

    for x in range(3,len(df9.columns)):
        value00 = (df9.iloc[0:count1,x] == '0/0').sum().sum()
        value11 = (df9.iloc[0:count1,x] == '1/1').sum().sum()

        if value00 >= value11 :
            df9.iloc[:,x].replace('0/0',0,inplace=True)
            df9.iloc[:,x].replace('1/1',2,inplace=True)
            df9.iloc[:,x].replace('0/1',1,inplace=True)
            df9.iloc[:,x].replace('1/0',1,inplace=True)
        else:
            df9.iloc[:,x].replace('1/1',0,inplace=True)
            df9.iloc[:,x].replace('0/0',2,inplace=True)
            df9.iloc[:,x].replace('0/1',1,inplace=True)
            df9.iloc[:,x].replace('1/0',1,inplace=True)

        avalue0 = (df9.iloc[0:count1,x] == 0).sum().sum()
        bvalue0 = (df9.iloc[count1:count1+count2,x] == 0).sum().sum()

        avalue1 = (df9.iloc[0:count1,x] == 1).sum().sum()
        bvalue1 = (df9.iloc[count1:count1+count2,x] == 1).sum().sum()

        caua = (2 * avalue0 + avalue1) / (2 * count1)
        caub = (2 * bvalue0 + bvalue1) / (2 * count2)

        freq_cutoff = float(freq_cutoff)
        if (caua-caub) >= freq_cutoff:
            resulttxt1 = df9.iloc[:, x:(x+1)]
            resulttxt = pd.concat([resulttxt, resulttxt1], axis=1)
        else:
            pass



    resulttxt.to_csv(out_file, sep='\t', index=False)
    os.system('rm test7.*')
    print('程序运行结束，请查看结果文件')



if __name__ == '__main__':
    main()













