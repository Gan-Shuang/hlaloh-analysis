import os
import sys
import time
import argparse
# import pymysql
from subprocess import call,getoutput,check_call
import pandas as pd
import numpy as np

# def purity_info(user,passwd,sample):
#     config = {
#           'host':'',
#           'port':'',
#           'user':'',
#           'password':'',
#           'database':'',
#           'charset':'utf8mb4',
#           'cursorclass':pymysql.cursors.Cursor
#           }
#     db = pymysql.connect(**config)
#     cursor = db.cursor()
#     cursor.execute("SELECT tumorPurity FROM # WHERE sampleId LIKE '%"+sample+"'")
#     purity_info = cursor.fetchone()
#     purity=purity_info[0]
#     return(purity)

def get_vaf_vcf(vcf):
    vaf_gene={}
    vcf_file=open(vcf,'r')
    for snv in vcf_file:
        snv=snv.rstrip()
        if not snv.startswith('#') and snv.split('\t')[6]=='PASS':
            gene=[]
            for ann in snv.split('\t')[7].split(';')[1].split(','):
                gene.append(ann.split('|')[3])
            gene=set(list(gene))
            for i in gene:
                if i=='TP53':
                    vaf=0
                    vaf_gene[str(gene)]=vaf
                    break
                else:
                    vaf=float(snv.split('\t')[7].split(';')[0].split('=')[1])
                    vaf_gene[str(gene)]=vaf
    if len(vaf_gene)==0:
        vaf_gene={'no_mutation':0}
    return(sorted(vaf_gene.values())[-1])

def run_hlaloh(normal_bam,out_dir,sample):
    mnt=str(out_dir).split('/')[1]
    work_path=sys.path[0]
    hlaloh_command=("singularity exec -B /"+mnt+":/"+mnt+" "+work_path+"/lohhla.sif "+
                    "Rscript "+work_path+"/lohhla/LOHHLAscript_setpath.R "+
                    "--HLAfastaLoc "+work_path+"/data/abc_complete.fasta "+
                    "--HLAexonLoc "+work_path+"/data/hla.dat "+
                    "--mappingStep TRUE "+
                    "--minCoverageFilter 10 "+
                    "--fishingStep TRUE "+"--cleanUp FALSE "+
                    "-k 31 "+
                    "--gatkDir "+work_path+"/picard/ "+
                    "--novoDir "+work_path+"/novocraft/ "+
                    "--patientId "+sample+" "+
                    "--normalBAMfile "+out_dir+"/tmp/link_bam/"+normal_bam+" "+
                    "--BAMDir "+out_dir+"/tmp/link_bam "+
                    "--outputDir "+out_dir+"/tmp "+
                    "--hlaPath "+out_dir+"/tmp/hlatype.txt "+
                    "--CopyNumLoc "+out_dir+"/tmp/solution.txt ")
    print(hlaloh_command)
    call(hlaloh_command,shell=True)

parser=argparse.ArgumentParser(description="Analyze HLA LOH")
parser.add_argument("-normal","-normalBam",dest="normal_bam_file",help="Input normal bam",required=True)
parser.add_argument("-tumor","-tumorBam",dest="tumor_bam_file",help="Input tumor bam",required=True)
parser.add_argument("-out","-outDir",dest="output_dir",help="Output dir",required=True)
parser.add_argument("-vcf","-vcf_file",dest="final_vcf_file",help="Sample somatic vcf(Get purity if no such information in database)")
parser.add_argument("-hla","-hlaType",dest="hlaType_file",help="HLA type information file",required=True)
args = parser.parse_args()
normal_bam=args.normal_bam_file
tumor_bam=args.tumor_bam_file
vcf=args.final_vcf_file
hla_input=args.hlaType_file
out_dir=args.output_dir

if __name__=="__main__":
    sample=tumor_bam.split('/')[-1].split('.')[0].split("-")[0]
    print(sample)
    call(["mkdir",out_dir])
    call(["mkdir",out_dir+"/tmp"])
    call(["mkdir",out_dir+"/tmp/link_bam"])
    # purity=purity_info('','',sample)
    purity=-99
    if purity<=0:
        if os.path.exists(vcf):
            purity=round(get_vaf_vcf(vcf),2)
    print('purity:{}'.format(purity))
    solution=open(out_dir+"/tmp/solution.txt","w")
    solution.write("Ploidy\ttumorPurity\ttumorPloidy\n"+tumor_bam.split('/')[-1].replace(".bam","")+"\t2\t"+str(purity)+"\t2\n")
    solution.close()
    hla_type=open(out_dir+"/tmp/hlatype.txt","w")
    ###hla format should be like 'hla_a|b|c_*_*'###
    for i in open(hla_input,"r"):
        i=i.rstrip()
        Allele=i.replace("HLA","hla").replace("A","a").replace("B","b").replace("C","c").replace("-","_").replace("*","_").replace(":","_")
        hla_type.write(Allele+'\n')
    hla_type.close()
    bam_ln_command=("ln -s "+tumor_bam+"* "+normal_bam+"* "+out_dir+"/tmp/link_bam")
    call(bam_ln_command,shell=True)
    ###run analysis###
    run_hlaloh(normal_bam.split("/")[-1],out_dir,sample)
    ###output to json(only HLA 1type)###
    HLA_A={}
    HLA_B={}
    HLA_C={}
    HLA_A["genesymbol"]="HLA-A"
    HLA_B["genesymbol"]="HLA-B"
    HLA_C["genesymbol"]="HLA-C"
    HLA_A["LOH"]="阴性"
    HLA_B["LOH"]="阴性"
    HLA_C["LOH"]="阴性"
    HLAlossPrediction=pd.read_csv(out_dir+"/tmp/"+sample+".10.DNA.HLAlossPrediction_CI.xls",sep="\t")
    ###P_value threshold set as 0.05###
    for index,data in HLAlossPrediction.iterrows():
        if HLAlossPrediction.loc[index]['HLA_A_type1'].startswith("hla_a"):
            if float(HLAlossPrediction.loc[index]['PVal_unique'])<=0.05:
                HLA_A["LOH"]="阳性"
        if HLAlossPrediction.loc[index]['HLA_A_type1'].startswith("hla_b"):
            if float(HLAlossPrediction.loc[index]['PVal_unique'])<=0.05:
                HLA_B["LOH"]="阳性"
        if HLAlossPrediction.loc[index]['HLA_A_type1'].startswith("hla_c"):
            if float(HLAlossPrediction.loc[index]['PVal_unique'])<=0.05:
                HLA_C["LOH"]="阳性"
    HLA_list=[HLA_A,HLA_B,HLA_C]
    print(HLA_list)
    result=open(out_dir+"/hlaloh.json","w")
    result.write(str(HLA_list)+'\n')
    result.close()
