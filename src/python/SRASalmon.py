import subprocess
import urllib.request
import xml.etree.ElementTree as ET
from sys import argv
import os
import json
import _thread

def run_thread(threadname, cmd):
    print("Started thread: ", threadname, ", Running: ", cmd)
    subprocess.run(' '.join(cmd), shell = True)

class accession:
    pid = ''
    paired = False
    runinfo = {}
    biosample = {}
    fastq = ''

    def __init__(self,pid):
        self.pid = pid

    def getRunInfo(self):
        url='https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&amp;db=sra&amp;rettype=runinfo&amp;term=' + self.pid
        u = urllib.request.urlopen(url)
        x = u.readlines()
        x0 = x[0].decode().strip().split(',')
        x1 = x[1].decode().strip().split(',')
        for i in range(0, len(x0)):
            self.runinfo[x0[i]] = x1[i]

    def getBioSample(self):
        if not self.runinfo:
            getRunInfo()
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=' + self.runinfo['BioSample']
        u = urllib.request.urlopen(url)
        tree = ET.parse(u)
        root = tree.getroot()
        for attr in root[0][5]:
            self.biosample[attr.attrib['attribute_name']] = attr.text
        bkeys = ['ecotype', 'tissue', 'dev_stage', 'sample_type', 'treatment', 'geo_loc_name',
                 'genotype', 'age']
        for k in bkeys:
            if not k in self.biosample:
                self.biosample[k]=''

    def getData(self):
        if not self.runinfo:
             getRunInfo()
        cmd = ['fasterq-dump', '-p']
        cmd.append('--split-files')
        cmd.append(self.pid)
        print(cmd)
        subprocess.run(cmd)
        self.fastq = self.pid + '.fastq'

    def salmon(self):
        infile = ''
        if os.path.exists(self.pid + '_2.fastq'):
            print("Mode: PE")
            cmd = ['singularity', 'exec',
                   '-B', '/mnt/picea/storage/reference/Drosophila-melanogaster/6.25/indices/salmon:/mnt/picea/storage/reference/Drosophila-melanogaster/6.25/indices/salmon:ro',
                   '/mnt/picea/projects/singularity/salmon.simg',
                   'salmon','quant','-p','8','-l','A','-1',self.pid + '_1.fastq',
                   '-2',self.pid + '_2.fastq','-i',
                   '/mnt/picea/storage/reference/Drosophila-melanogaster/6.25/indices/salmon/dmelr6.25.inx',
                   '--output', self.pid]
        else:
            if os.path.exists(self.pid + '_1.fastq'):
              infile = self.pid + '_1.fastq'
            else:
              infile = self.pid + '.fastq'
            print("Mode: SE")
            #cmd = ['srun','--mem', '16G', '-t','01-00:00:00','-c','8',
            cmd = ['singularity', 'exec',
                   '-B', '/mnt/picea/storage/reference/Drosophila-melanogaster/6.25/indices/salmon:/mnt/picea/storage/reference/Drosophila-melanogaster/6.25/indices/salmon:ro',
                   '/mnt/picea/projects/singularity/salmon.simg',
                   'salmon','quant','-l','A','-r',
                    infile,'-p','8','-i',
                   '/mnt/picea/storage/reference/Drosophila-melanogaster/6.25/indices/salmon/dmelr6.25.inx',
                   '--output', self.pid]
        print(cmd)
        subprocess.run(cmd)
        if os.path.exists(self.pid + '_2.fastq'):
          os.remove(self.pid + '_1.fastq')
          os.remove(self.pid + '_2.fastq')
        else:
           os.remove(infile)
        with open(self.pid + '.done', 'w') as done:
          done.write('1')

    def printTabular(self):
        r = self.runinfo
        b = self.biosample
        print(self.pid,
              r['Experiment'],r['Submission'],r['Sample'],r['BioSample'],r['BioProject'],r['LibraryLayout'],r['Platform'],
              r['spots'],r['avgLength'],r['bases'],r['InsertSize'],r['InsertDev'],r['spots_with_mates'],r['LibrarySource'],
              r['LibraryStrategy'],r['LibrarySelection'],r['TaxID'],r['size_MB'],r['ReleaseDate'],b['ecotype'],b['tissue'],
              b['dev_stage'],b['sample_type'],b['treatment'],b['geo_loc_name'],b['genotype'],b['age'],sep='\t')


if __name__ == '__main__':

    if os.path.exists(argv[1] + '.done'):
        print("Skipping", argv[1])
    else:
        print("Getting metadata...")
        x = accession(argv[1])
        try:
            x.getRunInfo()
            with open(argv[1] + "_runinfo.json", "w") as of:
                json.dump(obj=x.runinfo,fp=of)
            x.getBioSample()
            with open(argv[1] + "_biosample.json", "w") as of:
                json.dump(obj=x.biosample,fp=of)
        except:
            pass
        print("Getting data...")
        x.getData()
        print("Salmon...")
        x.salmon()
