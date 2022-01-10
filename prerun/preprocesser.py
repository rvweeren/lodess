#! /usr/bin/python3
import numpy as np
import multiprocessing as mp
import sys,os,glob
# import pyrap.tables as pt

'''
    This script downloads subbands from a staged html list.
    It automatically detects stalling and redownloads if needed
    Also, it splits the resulting files in separate folders

    =====================================================

    obviously, it would be better if the LTA would just work
    but that's to much to be asked
'''

class Pointing():
    def __init__(self,urls,lnum):
        self.urls = urls
        self.lnum = lnum

    def prerun(self):
        self.commands = []
        self.subs = []
        os.mkdir(self.lnum)
        for url in self.urls:
            self.commands.append(url[:-1])
            self.subs.append(url.split('/')[-1].split('_')[1])

    def check(self):
        downloaded = glob.glob(self.lnum+'/*')
        notexist = []
        for sub in self.subs:
            exists = False
            for gl in downloaded:
                if sub in gl:
                    exists = True
            if not exists:
                notexist.append(sub)
        self.commands = []
        if len(notexist) != 0:
            print(f'{self.lnum} still missing {len(notexist)} subbands...')
            submask = np.array([sb in notexist for sb in self.subs])
            urls_to_dl = self.urls[submask]
            for url in urls_to_dl:
                self.commands.append(url[:-1])
        else:
            print(f'{self.lnum} complete!')

    def renameuntar(self):
        '''
            Does both renaming/untarring/deleting
            But also checks if there are invalid tars, which have to be redownloaded again
        '''
        os.chdir(self.lnum)
        baddl = []
        for filename in glob.glob("*SB*.tar"):
            outname=filename.split("%")[-1]
            os.rename(filename, outname)
            exitval = os.system('tar -xvf '+outname)
            if exitval == 0:
                os.system('rm -r '+outname )
            else:
                baddl.append(filename.split('_')[-3])
                os.rename(outname, filename)
        os.chdir('..')
        self.commands = []
        if len(baddl) != 0:
            print(f'{self.lnum} still missing {len(baddl)} subbands...')
            submask = np.array([sb in baddl for sb in self.subs])
            urls_to_dl = self.urls[submask]
            for url in urls_to_dl:
                self.commands.append(url[:-1])
            print(self.commands)
        else:
            print(f'Untarring {self.lnum} complete!')

    def checkrun(self):
        '''
            Downloads iteratively, and then checks if all files are ok
        '''
        while len(self.commands) != 0:
            os.chdir(self.lnum)
            for cmd in self.commands:
                command = f'wget -c --timeout=10 --tries=999 {cmd}'
                os.system(command)
            os.chdir('..')
            self.check()
            if len(self.commands) == 0:
                self.renameuntar()


def process_pointing(run):
    run.prerun()
    run.checkrun()

def process_html(file):
    data = []
    with open(file,'r') as handle:
        for line in handle:
            data.append(line)
    lnums = np.array([lin.split('/')[-1].split('_')[0] for lin in data])
    threads = []
    for lnum in np.unique(lnums):
        matches = np.array(data)[lnums == lnum]
        run = Pointing(matches,lnum)
        thread = mp.Process(target=process_pointing,args=(run,),daemon=True)
        thread.start()
        threads.append(thread)
    for i in threads:
        i.join()

if __name__  == '__main__':
    process_html(sys.argv[1])
