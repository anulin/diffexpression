import glob
import os
import codecs

humanfactors={}
mousefactors={}
mtissues={}
htissues={}
factors=set()
HTOM={}
n=0
with open("HOM_MouseHumanSequence.rpt.txt") as file:#соответствие генов человека и мыши
    for line in file:
        if line.split()[1]=="mouse":
            mouse=line.split()[3]
            n=line.split()[0]
        if line.split()[1] == "human" and line.split()[0]==n:
            HTOM[line.split()[3]]=mouse
with open("факторы с id.txt") as file:
    for line in file:
        if(line.split()[0] in HTOM):
            factors.add(line.split()[0])
            factors.add(HTOM[line.split()[0]])
GoodTissues={}
n=0
name=""

humangenes={}

with open("белки.txt") as file:#hgnc id в symbol и symbol в ensg
    for line in file:
        if line.split()[1] in HTOM:
            humangenes[line.split()[0]] = line.split()[1]
mousegenes={}
with open("M20.txt") as file:# uniprot B symbol
    for line in file:
        if(line.split()[1]in HTOM.values()):
            mousegenes[line.split()[0]] = line.split()[1]
hcis={}#конветрация нормальных названий факторов транскрипции
mcis={}
with open("Human cistromF.txt") as file:
    for line in file:
        if(line.split()[1]in HTOM):
            hcis[line.split()[0]] = line.split()[1]
for filename in hcis:
    if ("cistrome2genes_hg19\\hg19.gencode_v19." + filename + '.tsv') in glob.glob(
            os.path.join("cistrome2genes_hg19", '*.tsv')):
        with open("cistrome2genes_hg19/hg19.gencode_v19."+filename+'.tsv') as file:
            humanfactors[hcis[filename]]=set()

with open("mouse cistromF.txt") as file:
    for line in file:
        if(line.split()[1]in HTOM.values()):
            mcis[line.split()[0]] = line.split()[1]

for filename in mcis :
    if ("cistrome2genes_mm9\\mm9.gencode_vM1." + filename + '.tsv') in glob.glob(
            os.path.join("cistrome2genes_mm9", '*.tsv')):
        with open("cistrome2genes_mm9/mm9.gencode_vM1."+filename+'.tsv') as file:
            mousefactors[mcis[filename]]=[]

exprf={}

positions=[]
names=[]
hgcids=set()
with open("gene_expressions.tsv") as file:
    i=0
    for line in file:
        while i< len(line.split()):
            if len(line.split()[i].split("."))>3 :# ДОбавить фильтр больных тканей!!!!!!
                positions.append(i)
            i+=1

        if line.split()[3] in humangenes and humangenes[line.split()[3]] in humanfactors:
            if(humangenes[line.split()[3]] not in exprf):
                exprf[humangenes[line.split()[3]]]=0
            for j in positions:
                if float(line.split()[j])>exprf[humangenes[line.split()[3]]]:
                    exprf[humangenes[line.split()[3]]]=float(line.split()[j])

positionsm=[]
namesm=[]
exprmf={}

with open("gene_expressionsM.tsv") as file:
    i = 0
    for line in file:
        while i < len(line.split()):
            if len(line.split()[i].split("."))>3:#
                positionsm.append(i)
            i += 1

        if line.split()[3] in mousegenes and mousegenes[line.split()[3]] in mousefactors:
            if(mousegenes[line.split()[3]] not in exprmf):
                exprmf[mousegenes[line.split()[3]]] = 0
            for j in positionsm:
                if float(line.split()[j]) > exprmf[mousegenes[line.split()[3]]]:
                    exprmf[mousegenes[line.split()[3]]]=float(line.split()[j])

mth={}
i=0
print(len(exprf), len(exprmf))
print(exprmf)
with open("MaxFactExprs.txt",'w') as file:
    for fact in exprmf:
        file.write(fact+' '+str(exprmf[fact])+'\n')
    for fact in exprf:
        file.write(fact + ' ' +str( exprf[fact]) + '\n')
