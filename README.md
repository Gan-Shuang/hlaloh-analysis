# hlaloh-analysis
HLA-LOH detect method for most NGS data.  
## Introduction
hlaloh-analysis based on https://github.com/mskcc/lohhla. suitable for most NGS data covered HLA area.
## Install
### build singularity evironment 
Reqired R tools and samtools,bedtools in this env 
```
cd hlaloh-analysis  
singularity build lohhla.sif docker://ganshuang0925/samtools_bedtools  
```
### install jellyfish
package was offered
```  
cd hlaloh-analysis  
tar -zxvf jellyfish-1.1.10.tar.gz  
cd jellyfish-1.1.10  
./configure --prefix=${PWD}/../jellyfish  
make  
make install  
```
  
### change work path
change R script tools path
```
cd hlaloh-analysis  
sed "s#MAIN#${PWD}#g" ${PWD}/lohhla/LOHHLAscript.R > ${PWD}/lohhla/LOHHLAscript_setpath.R  
chmod a+x ${PWD}/lohhla/LOHHLAscript_setpath.R  
```
### unzip reference
```
cd hlaloh-analysis  
cd data  
gizp -d *gz
```
## Method
```
python3 hlaloh.py -h
optional arguments:
  -h, --help            show this help message and exit
  -normal NORMAL_BAM_FILE, -normalBam NORMAL_BAM_FILE
                        Input normal bam
  -tumor TUMOR_BAM_FILE, -tumorBam TUMOR_BAM_FILE
                        Input tumor bam
  -out OUTPUT_DIR, -outDir OUTPUT_DIR
                        Output dir
  -vcf FINAL_VCF_FILE, -vcf_file FINAL_VCF_FILE
                        Sample somatic vcf(Get purity if no such information in database)
  -hla HLATYPE_FILE, -hlaType HLATYPE_FILE
                        HLA type information file
```
note：purity information was counted from VAF from VCF file（or loaded by modificate py script）  
## Test
All inputs should be absolute path  
```
cd hlaloh-analysis 
python3 ${PWD}/hlaloh.py -tumor ${PWD}/test_example/test-tumor.bam -normal ${PWD}/test_example/test-normal.bam -hla ${PWD}/test_example/test-hlatype -vcf ${PWD}/test_example/test.vcf -out ${PWD}/test_out
```
## Result
> hlaloh.json
```
[{'genesymbol': 'HLA-A', 'LOH': '阴性'}, {'genesymbol': 'HLA-B', 'LOH': '阳性'}, {'genesymbol': 'HLA-C', 'LOH': '阴性'}]
```
> tmp dir
>> test.10.DNA.HLAlossPrediction_CI.xls
