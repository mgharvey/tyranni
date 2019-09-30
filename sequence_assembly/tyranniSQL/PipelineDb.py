#!/usr/bin/env python2
##useful

from __future__ import print_function
import itertools, os, io
#import datetime
from tyranniSQL.helpers import StdoutWordWrap,MessageTicker,TickerWriter
from tyranniSQL.Stage3LociFile import Stage3LociFile
from tyranniSQL.stage4 import Stage4SamplesFile
from tyranniSQL.stage5 import Stage5SamplesFile
from Bio import SeqIO
#from numpy import mean
import mysql.connector
import re,os.path
#import sys,string
import platform
#import openpyxl,lxml
#from openpyxl import load_workbook
#from mysql.connector.errors import IntegrityError
from collections import namedtuple
from itertools import imap
import operator, gzip

"""

Name: PipelineDb.py

Author: Graham E Derryberry
Date: 11 September 2016

Description:
Tools for connecting to and interacting with the database.

"""

#C4CMUACXX_s1_1_RapidGenomics_10_SL63193.fastq.gz
#C6JLDANXX_s5_1_RapidGenomics_9-TruseqHTDual_D503_SL113032.fastq.gz
rawRE=re.compile('((..).(..)..XX)_s(.)_(.)_RapidGenomics_([0-9]+).*?_SL([0-9]+)\\.fastq\\.(.+)')
codeRE=re.compile('((..).(..)..XX)_s(.)_(.)_RapidGenomics_([0-9]+).*?_SL([0-9]+)\\.fastq')
tmpSIDre=re.compile('[^0-9]+')

JobKey=namedtuple("JobKey",["Batch","name"])

FilepathRow=namedtuple("FilepathRow",["inode","Batch","job_name","file_name","volume","abspath"])
STAGE1=u'stage1'
reS1=re.compile('^stage1$',flags=re.I)
reS1Batch=re.compile('^stage1|^uce',flags=re.I)
def is_stage1(l):
    return(reS1.match(l) is not None)
STAGE2=u'stage2'
reS2=re.compile('^stage2$',flags=re.I)
reS2Batch=re.compile('^stage2',flags=re.I)
def is_stage2(l):
    return(reS2.match(l) is not None)
STAGE3=u'stage3'
reS3=re.compile('^stage3$',flags=re.I)
reS3Batch=re.compile('^stage3|^uce',flags=re.I)
def is_stage3(l):
    return(reS3.match(l) is not None)
STAGE4=u'stage4'
reS4=re.compile('^stage4$',flags=re.I)
reS4Batch=re.compile('^stage4|^uce',flags=re.I)
def is_stage4(l):
    return(reS4.match(l) is not None)
STAGE5=u'stage5'
reS5=re.compile('^stage5$',flags=re.I)
reS5Batch=re.compile('^stage5|^uce',flags=re.I)
def is_stage5(l):
    return(reS5.match(l) is not None)

def get_rgID(mG):
    vals=list(mG)
    vals[5]=int(vals[5])
    return "rg{2}{3}{5:02d}_{1}_T".format(*vals)


def str_rVal(rVal):
    if rVal is None:
        return(None)
    if type(rVal) in [str,unicode]:
        return(rVal)
    return(str(rVal))

class PipelineDb:
    def __init__(self,conn=None,cxn_args=None,dry_run=False,rollback=False,clone=False,ticker=None):
        self.dry_run=dry_run
        self.conn=conn
        self.cxn_args=cxn_args
        self.rollback=rollback
        if clone:
            return
        if ticker is None:
            self.ww=StdoutWordWrap()
        else:
            self.ww=ticker
        self.ticker=MessageTicker(self.ww)
        self.ticker_out=TickerWriter(self.ticker)
        self.locus_ref=None
        self.loci=None
        self.location=detect_location()
        if conn is not None:
            print("Pipeline connection made to {}.".format(self.conn.database),file=self.ticker)
    def clone_conn(self):
        alt=PipelineDb(cxn_args=self.cxn_args,dry_run=self.dry_run,rollback=self.rollback,clone=True)
        alt.ww=self.ww
        alt.ticker=MessageTicker(self.ww)
        alt.locus_ref=self.locus_ref
        alt.loci=self.loci
        alt.location=self.location
        alt.ticker_out=TickerWriter(alt.ticker)
    def open_connection(self):
        self.conn=mysql.connection.connect(**self.cxn_args)
        print("Pipeline connection made to {}.".format(self.conn.database),file=self.ticker)
    def desc_table(self,t):
        res=zip(*(self.conn.cmd_query('''Select * from {0} limit 1'''.format(t))['columns']))[0]
        self.conn.get_rows()
        return(res)
    def close(self):
        print("Closing pipeline connection.",file=self.ticker)
        self.conn.close()
        self.conn=None
    def named_cursor(self):
        if self.conn is None:
            self.open_connection()
        return self.conn.cursor(named_tuple=True)
    def parse_loci_file(self,filepath):
        "Encapsulate a Stage 3 file"
        return Stage3LociFile(filepath,self)
    def select_locus_sequences(self,locus):
        "Prepare to download Stage 4 file for locus"
        if self.locus_ref is None:
            return Stage4SamplesFile(self,locus_name=locus)
        return Stage4SamplesFile(self,locus_name=locus,locus_uid=self.locus_ref[locus])
    def select_locus_alignments(self,locus,group):
        "Prepare to download Stage 4 file for locus"
        if self.locus_ref is None:
            return Stage5SamplesFile(self,locus_name=locus,group_name=group)
        return Stage5SamplesFile(self,locus_name=locus,locus_uid=self.locus_ref[locus],group_name=group)
        #raise NotImplementedError("Still no stage 4 code.")
    def parse_locus_alignments(self,filepath):
        "Interact with locus alignment in database"
        #raise NotImplementedError("Still no stage 5 code.")
        #if self.locus_ref is None:
        #    return Stage5SamplesFile(self,locus_name=locus)
        #return Stage5SamplesFile(self,locus_name=locus,locus_uid=self.locus_ref[locus])
        return Stage5SamplesFile(self,fp=filepath)
        #raise NotImplementedError("Still no stage 4 code.")
    def load_locus_dict(self):
        "Downloads and stores the mapping of loci names to ids."
        c=self.conn.cursor()
        c.execute('''select locus_name,locus_uid from loci order by locus_uid''')
        res=c.fetchall()
        self.locus_ref=dict(res)
        self.loci=self.locus_ref.keys()
        self.locus_dict=dict(zip(self.locus_ref.values(), self.locus_ref.keys()))
    def get_sample_dict(self,stage='stage4'):
        if stage is 'stage4':
            c=self.conn.cursor()
            c.execute('''select sample_uid,SampleName 
            from sample_aliases where  prefer_align=b'1' order by sample_uid;''')
            res=c.fetchall()
            return(dict(res))
        raise NotImplementedError("No sample name dictionary yet for stage",stage)
    def get_sample_rev_dict(self,stage='stage4'):
        if stage is 'stage4':
            c=self.conn.cursor()
            c.execute('''select SampleName,sample_uid
            from sample_aliases where  prefer_align=b'1' order by sample_uid;''')
            res=c.fetchall()
            return(dict(res))
        raise NotImplementedError("No sample name dictionary yet for stage",stage)
    def get_locus_uid(self,locus_name):
        "Simple method to fetch id number from database.  Providing an iterable will produce a dictionary."
        if self.locus_ref is None:
            self.load_locus_dict()
        return(self.locus_ref[locus_name])
    def get_locus(self,locus_uid):
        "Simple method to fetch id number from database.  Providing an iterable will produce a dictionary."
        if self.locus_ref is None:
            self.load_locus_dict()
        return(self.locus_dict[locus_uid])
    def load_job_file_from_cache(self,filepath):
        "Use cached filepath to instantiate a JobFile"
        fp=os.path.abspath(filepath)
        c=self.named_cursor()
        c.execute("""SELECT dir_path, stage, batch from batch_path_cache 
            where %s regexp binary dir_path and stage is not null and batch is not null 
            and system = %s""",
            (os.path.abspath(filepath),self.location))
        res=c.fetchall()
        if c.rowcount is 0:
            raise Exception("Uncached batch for",filepath)
        if c.rowcount > 1:
            print("Warning:",filepath,"matched to multiple batches in file cache")
        for b_row in res:
            b=Batch(self,b_row.dir_path)
            if re.match("stage[123]",b.stage):
                f_key=break_path(fp,root=b.base_path)
            else:
                f_key=seq_break_path(fp,root=b.base_path)
            j=b.job_file(f_key)
            j.pull()
            if j.file_uid is not None:
                return j
        #raise Exception("Uncached file",filepath)
        print("Warning, uncached file",filepath)
        j.push()
        return j
    def load_file_uid_from_cache(self,filepath,force=False):
        "Use cached filepath to retrieve file_uid"
        fp=os.path.abspath(filepath)
        c=self.named_cursor()
        c.execute("""SELECT file_uid from file_cache 
            where %s = abspath and system = %s""",
            (fp,self.location))
        res=c.fetchall()
        if c.rowcount is 0:
            if not force:
                return None
            jF=self.load_job_file_from_cache(filepath)
            jF.pull();jF.push()
            return jF.file_uid
        if c.rowcount > 1:
            print("Warning:",filepath,"matched to multiple files in file cache")
        return res[0].file_uid
    def tick(self,text='X'):
        print(text,end='',file=self.ticker_out)
    def flush_ticks(self):
        self.ticker.flush()
    def push_ticks(self):
        self.ticker.push()

def detect_location():
    "Figure out whether we are running on cypress, eeb-globus, or someplace else."
    if re.match('^cypress', platform.node()) is not None:
        return('cypress')
    if re.match('^eeb-globus', platform.node()) is not None:
        return('eeb-blobus')
    return(platform.node())

def get_inode_tuple(filepath):
    "Return the tuple (inode, filesystem, filepath)"
    if not os.path.exists(filepath):
        raise Exception("{} does not exist!".format(filepath))
    if os.path.isdir(filepath):
        raise Exception("{} is a directory!".format(filepath))
    return((os.stat(filepath).st_ino,get_path_volume(filepath),os.path.abspath(filepath)))


if platform.system() in ['Linux',u'Linux','Darwin',u'Darwin']:
    reM=re.compile('.*? on (/.*?)( type [a-zA-Z0-9]+) \([^)]+\)')
    vols=[ m.groups()[0] for m in map(lambda TTT: reM.match(TTT),  os.popen('mount -v') ) if m is not None]
    vols.sort(reverse=True)
    reVols=re.compile('^(' + '|'.join(vols)+')')
else:
    reVols=re.compile('^(.):')

def get_path_volume(path):
    "Returns the volume on which the target of path resides"
    return(reVols.match(os.path.realpath(path)).groups()[0])

def break_path(filepath,root=None,batch=None,job=None):
    "Breaks a path up into key parts.  root is the base of the batch directory."
    fp=filepath
    filepath=os.path.normpath(filepath)
    if root is not None:
        filepath = filepath.replace(os.path.join(root,''),'')
    parts=filepath.split(os.sep)
    if job is not None and batch is None:
        raise Exception("Setting the job name with out the batch name is not supported")
    if job is not None:
        parts = [job] +parts
    if batch is not None:
        parts = [batch] +parts
    iNd, fs, fullp= get_inode_tuple(fp)
    return( FilepathRow(iNd,unicode(parts[0]),unicode(parts[1]),unicode(os.path.join(*parts[2:])), unicode(fs),unicode(fullp)))

def seq_break_path(filepath,root=None,batch=None,job=None):
    "Breaks a path up into key parts.  root is the base of the batch directory."
    fp=filepath
    filepath=os.path.normpath(filepath)
    if root is not None:
        filepath = filepath.replace(os.path.join(root,''),'')
    parts=filepath.split(os.sep)
    if job is not None and batch is None:
        raise Exception("Setting the job name with out the batch name is not supported")
    if job is None:
        tmp = reSEQfile.search(os.path.basename(filepath));
        if tmp is None: 
            raise Exception("Not a uce sequence file!")
        job = tmp.groups()[0]
    if batch is not None:
        parts = batch.split(os.sep) +parts
    iNd, fs, fullp= get_inode_tuple(fp)
    return( FilepathRow(iNd,unicode(os.path.join(*parts[0:2])),unicode(job),unicode(os.path.join(*parts[2:])), unicode(fs),unicode(fullp)))

def db2jk(row):
    return(JobKey(Batch=unicode(row.Batch),name=unicode(row.job_name)))

def db2jkDict(rows):
    r1,r2=itertools.tee(rows)
    return(dict(itertools.izip(itertools.imap(db2jk,r1),itertools.imap(lambda r: r.job_uid,r2))))

## files of interest
#examples
#     rgUM667_C4_T-READ1.fastq.gz
#     stage2-manakin/Antilo_gal_L14634/Antilo_gal_L14634.fasta
#     stage3-manakin/Antilo_gal_L14634/Antilo_gal_L14634_contigs.fasta
#     stage3-manakin/Antilo_gal_L14634/Antilo_gal_L14634_loci.fasta
#     C4CMUACXX_s1_1_RapidGenomics_20_SL63200.fastq.gz
reI=re.compile('rg....._.._.-READ..fastq.gz|RapidGenomics|[A-Z][^.]*.fasta$|contigs.fa$|cns.fq.gz$|full_depth.gz$|depth.txt$')
reSEQfile=re.compile('(uce-[^.]*|[^_.]*_exon[^.]*).fasta$');

def file_of_interest(filename):
    "Is this file worth processing?"
    filename=os.path.basename(filename)
    return(reI.search(filename) is not None)

def seq_file_of_interest(filename):
    "Is this file worth processing?"
    filename=os.path.basename(filename)
    return(reSEQfile.search(filename) is not None)


def walk_batch(batch_name,func=print):
    "Determine interesting files and pass key parts to func"
    for p,d,f in os.walk(batch_name): 
        f=[i for i in f if file_of_interest(i) ]
        if len(f) is 0:
            continue
        #print( f)
        fp=map(lambda TTT: os.path.join(p,TTT),f)
        #print(fp)
        for i in map(break_path,fp):
            func(i)

def walk_seq_batch(batch_name,func=print):
    "Determine interesting files and pass key parts to func"
    for p,d,f in os.walk(batch_name): 
        f=[i for i in f if seq_file_of_interest(i) ]
        if len(f) is 0:
            continue
        #print( f)
        fp=map(lambda TTT: os.path.join(p,TTT),f)
        #print(fp)
        for i in map(seq_break_path,fp):
            func(i)
# 
# 
# rows=list()
# walk_batch('uce-manakin',func=lambda TTT: rows.append(TTT))
# print(*rows[0:10],sep='\n')
# 
# rows=list()
# walk_batch('stage2-manakin',func=lambda TTT: rows.append(TTT))
# print(*rows[0:10],sep='\n')
# 
# rows=list()
# walk_batch('stage3-manakin',func=lambda TTT: rows.append(TTT))
# print(*rows[0:10],sep='\n')
# 
# 

## sql blips
all_jobs=''' select distinct Batch,job_name from jobs'''
batch_jobs=''' select distinct Batch,job_name,job_uid from jobs where batch = %s'''
group_jobs=''' select distinct Batch,job_name,job_uid from jobs where batchgroup = %s'''
find_batchpath='''select * from batch_path_cache 
where dir_path=%s and system=%s'''
find_batchjob='''select * from jobs 
where Batch=%s limit 1'''
find_batchstat='''select * from (
SELECT distinct Batch,BatchGroup,Branch,'stage2' as stage FROM stage2stats union
SELECT distinct Batch,BatchGroup,Branch,'stage3' as stage FROM stage3stats
) as sX where Batch=%s limit 1'''



def unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in itertools.ifilterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

def unique_justseen(iterable, key=None):
    "List unique elements, preserving order. Remember only the element just seen."
    # unique_justseen('AAAABBBCCDAABBB') --> A B C D A B
    # unique_justseen('ABBCcAD', str.lower) --> A B C A D
    return imap(next, imap(operator.itemgetter(1), itertools.groupby(iterable, key)))


class PipelineObj(object):
    "Represents an item with a unique id and a PipelineDb Connection"
    def __init__(self, pldb,auto_commit=False):
        if issubclass(type(pldb),PipelineObj):
            self.pldb=pldb.pldb
        else:
            self.pldb=pldb
        self.auto_commit=auto_commit
    def check_connection(self):
        if self.pldb.conn is None:
            self.pldb.open_connection()
    def pull(self, auto_open=True):
        "Stub for objects which have information cached in the db. Opens the connection."
        if auto_open:
            self.check_connection()
    def push(self):
        "Automatically commits changes if auto_commit is True and the db will commit."
        if not self.auto_commit: return;
        if self.pldb.dry_run: return;
        if self.pldb.conn is None:
            print("Warning, there is no connection to the database.")
            return;
        if self.pldb.rollback: return;
        self.pldb.conn.commit();
    def named_cursor(self):
        self.check_connection()
        return self.pldb.conn.cursor(named_tuple=True)
    def prepared_cursor(self):
        self.check_connection()
        return self.pldb.conn.cursor(prepared=True)
    def tick(self,text='#'):
        self.pldb.tick(text)
    def flush_ticks(self):
        self.pldb.flush_ticks()
    def push_ticks(self):
        self.pldb.flush_ticks()

## not sure if I should integrate or not
class Batch(PipelineObj):
    "Represents a batch of jobs that may or may not have been added to the database"
    def __init__(self, pldb, batchpath, batchgroup=None,batchbranch=None):
        PipelineObj.__init__(self,pldb)
        self.path=batchpath
        self.group=batchgroup
        self.batch_uid=None
        self.stage=None
        self.branch=batchbranch
        self.name=os.path.basename(batchpath)
        self.base_path=os.path.dirname(os.path.abspath(batchpath))
        if re.match('exon$|uce[0-9]$',self.name):
            self.name=os.path.join(os.path.basename(os.path.dirname(batchpath)),
                os.path.basename(batchpath))
            self.base_path=os.path.dirname(self.base_path)
        self.job_ref=None
        self.load_info()
    def load_info(self):
        "Figure out if batch is known"
        #Simple Test: is batchpath in cache?
        c=self.named_cursor()
        #print(find_batchpath,(os.path.realpath(self.path),self.pldb.location))
        c.execute(find_batchpath,(os.path.realpath(self.path),self.pldb.location))
        rows=c.fetchall()
        if c.rowcount >0 :
            row=rows[0]
            print("Found cached entry:",row)
            if re.match(row.Batch,self.name) is None:
                raise Exception("Integrity Error!")
            if self.group is not None and re.match(row.Batchgroup,self.group) is None:
                print("Warning: {} != {}".format(row.Batchgroup,self.group))
            if self.branch is not None and re.match(row.Branch,self.branch) is None:
                raise Exception("Integrity Error!")
            self.group=row.Batchgroup
            self.batch_uid=row.batch_uid
            self.stage=row.stage
            self.branch=row.Branch
            self.load_db_jobs()
            print('=',end='',file=self.pldb.ww)
            return
        #Check records: is batch mentioned in a known job?
        c.execute(find_batchjob,(self.name,))
        rows=c.fetchall()
        if c.rowcount >0 :
            row=rows[0]
            print("Found job entry:",row)
            if re.match(row.Batch,self.name) is None:
                raise Exception("Integrity Error!")
            if self.group is not None and re.match(row.BatchGroup,self.group) is None:
                print("Warning: {} != {}".format(row.BatchGroup,self.group))
            if self.branch is not None and re.match(row.Branch,self.branch) is None:
                raise Exception("Integrity Error!")
            self.group=row.BatchGroup
            self.stage=row.stage
            self.branch=row.Branch
            self.load_db_jobs()
        #Check records: is batch mentioned in a known stat file?
        c.execute(find_batchstat,(self.name,))
        rows=c.fetchall()
        if c.rowcount >0 :
            row=rows[0]
            print("Found stat entry:",row)
            if re.match(row.Batch,self.name) is None:
                raise Exception("Integrity Error!")
            if self.group is not None and re.match(row.BatchGroup,self.group) is None:
                print("Warning: {} != {}".format(row.BatchGroup,self.group))
            if self.branch is not None and re.match(row.Branch,self.branch) is None:
                raise Exception("Integrity Error!")
            self.group=row.BatchGroup
            self.stage=row.stage
            self.branch=row.Branch
        if self.group is None:
            self.group=os.path.basename(os.path.dirname(os.path.abspath(self.path)))
        if self.stage is None:
            m=re.match('^(stage[^-]*)-',self.name)
            if m is not None:
                self.stage=unicode(m.groups()[0])
            else:
                m=re.match('^uce-',self.name)
                if m is not None:
                    self.stage=u'stage1'
                else:
                    raise Exception("Could not determine stage of {}.".format(self.path))
        if self.branch is None:
            self.branch=u'normal'
        if self.job_ref is None:
            self.job_ref=dict()
        ins='''INSERT INTO batch_path_cache (system,batchgroup,batch,stage,branch,dir_path,base_path)
        values ({})'''.format(','.join(['%s']*7))
        vals=(self.pldb.location,self.group,self.name,self.stage,self.branch,os.path.abspath(self.path),self.base_path)
        if self.pldb.dry_run:
            #print(ins)
            print(vals[0:5])
            print('N',end='',file=self.pldb.ww)
        else:
            c.execute(ins,vals)
            self.batch_uid=c.lastrowid
            print('+',end='',file=self.pldb.ww)
#     def reload_keys(self):
#         c=self.pldb.conn.cursor(named_tuple=True)
    def iter_files(self):
        "Generator to iterate over all iteresting files and return paths"
        if re.match('stage[123]',self.stage ) is not None:
            for p,d,f in os.walk(self.path): 
                f=[i for i in f if file_of_interest(i)]
                d0=[i for i in d  if i.startswith('auto_data')]
                for i in d0:
                    d.remove(i)
                if len(f) is 0:
                    continue
                # if os.path.basename(p).startswith('auto_data'):
                    # continue
                #print( f)
                fp=map(lambda TTT: os.path.join(p.replace(self.path,self.name) ,TTT),f)
                #print(fp)
                for i in map(break_path,fp):
                    yield i
        elif re.match('stage[45]',self.stage ) is not None:
            for p,d,f in os.walk(self.path): 
                fS=(i for i in f if seq_file_of_interest(i))
                fp=itertools.imap(lambda TTT: os.path.join(p.replace(self.path,self.name) ,TTT),fS)
                #print(fp)
                for i in itertools.imap(seq_break_path,fp):
                    yield i
        else:
            raise Exception("Unknown stage: ",self.stage)
    def iter_jobs(self):
        "Generator to iterate over all jobs with iteresting files and return unique keys"
        return (unique_justseen(imap(db2jk,self.iter_files())) )
    def iter_batch_jobs(self):
        "Generator to iterate over all jobs with iteresting files and return unique objects"
        return(imap (lambda TTT: BatchJob(self,TTT),self.iter_jobs()))
    def iter_job_files(self):
        "Generator to iterate over all iteresting files and return JobFile objects"
        return (imap (lambda TTT: JobFile(self,TTT),self.iter_files()))
    def load_db_jobs(self):
        "Load all previously assigned job_uid back references"
        c=self.pldb.conn.cursor(named_tuple=True)
        c.execute(batch_jobs,(self.name,))
        self.job_ref=db2jkDict(c.fetchall())
    def register_job(self,jk,job_status=None):
        "Insert new job row into db"
        c=self.pldb.conn.cursor(named_tuple=True)
        if self.pldb.dry_run:
            print(self.insJob)
            print((jk.job_name,job_status))
            print('N',end='',file=self.pldb.ww)
        else:
            c.execute(self.insJob,(jk.job_name,job_status))
            self.job_ref.put(jk,c.lastrowid)
            print('+',end='',file=self.pldb.ww)
    def job_file(self, filekey):
        return JobFile(self,filekey)
    def tick(self,text='%'):
        self.pldb.tick(text)
    def push_ticks(self):
        self.pldb.push_ticks()


# def job_status():
#     """0: Not yet run
# 1: Succeeded
# -1: Failed complete
# -2: Results rejected
# 2: Ghost job (output only)
# 3: archived data (output and stats)
# 4: archived log (stats only)
# 5: archived reference (Job must have occurred but no other information exists)
# 6: Soft link (use ghost_ref to find and link to files)"""
#     pass

# def job_file_purpose():
#     """NULL: not yet determined
# 0: INPUT: used by job to run, not used in a future stage
# 1: OUTPUT: created by job for use in a future stage
# 2: CARRY: used by job to run, and in future stages
# 3: LOG: intended to carry extra information not yet stored in database
# 4: EXTRA: created by job but not used by pipeline"""
#     pass

getJobByUID='''SELECT job_uid,ghost_ref,job_name,stage,Branch,Batch,BatchGroup,job_status
FROM jobs where job_uid=%s'''
getJobByRef='''SELECT job_uid,ghost_ref,job_name,stage,Branch,Batch,BatchGroup,job_status
FROM jobs where job_name=%(name)s and batch=%(Batch)s '''
insJobRef='''INSERT INTO jobs (ghost_ref,job_name,stage,Branch,Batch,BatchGroup,job_status)
VALUES ({})'''.format(','.join(['%s']*7))

getSbyUID={u'stage1':'''SELECT * FROM stage1stats where job_uid=%s''',
u'stage2':'''SELECT * FROM stage2stats where job_uid=%s''',
u'stage3':'''SELECT * FROM stage3stats where job_uid=%s'''}
getSbyRef={u'stage1':'''SELECT * FROM stage1stats join jobs using (job_uid) where job_name=%(name)s and batch=%(Batch)s ''',
u'stage2':'''SELECT * FROM stage2stats where SampleName=%(name)s and batch=%(Batch)s ''',
u'stage3':'''SELECT * FROM stage3stats where SampleName=%(name)s and batch=%(Batch)s '''}
updateSJUID={
STAGE2: '''UPDATE stage2stats set job_uid=%(UID)s where SampleName=%(name)s and batch=%(Batch)s''',
STAGE3:  '''UPDATE stage3stats set job_uid=%(UID)s where SampleName=%(name)s and batch=%(Batch)s'''}

class BatchJob(PipelineObj):
    "Simple class to wrap around a job directory"
    def __init__(self, batch, jobkey):
        "Load up parts"
        PipelineObj.__init__(self,batch)
        self.batch   = batch   #Batch object
        self.jobkey  = jobkey  #JobKey named tuple
        self.job_uid = None
        self.ghost_ref = None
        self.job_status = None
        self.modified = True
        if jobkey in batch.job_ref.keys():
            self.job_uid=batch.job_ref[jobkey]
#        self.cache_uid = None
    def pull(self):
        "Get row, if it exists, from db"
        c = self.batch.pldb.named_cursor()
        if self.job_uid is not None:
            c.execute(getJobByUID,(self.job_uid,))
        else:
            c.execute(getJobByRef,self.jobkey.__dict__)
        rows=c.fetchall()
        # print(*rows,sep='\n')
        if c.rowcount <1:
            return
        if c.rowcount >1:
            print("Warning, multiple job match, using just the first")
       	self.modified =	False
        if self.job_uid is None and rows[0].job_uid is not None:
            self.job_uid=rows[0].job_uid
            self.batch.job_ref[self.jobkey]=self.job_uid
        if self.job_uid is None:
            self.modified = True
        if self.job_status is None and rows[0].job_status is not None:
            if rows[0].job_status <0 or rows[0].job_status>1:
                self.job_status=rows[0].job_status
    def pull_stats(self):
        "Get row, if it exists, from db"
        if self.batch.stage not in getSbyRef.keys():
            return
        c = self.batch.pldb.named_cursor()
        #First by ref
        c.execute(getSbyRef[self.batch.stage],self.jobkey.__dict__)
        rows=c.fetchall()
        # print(*rows,sep='\n')
        if self.job_uid is not None and c.rowcount<1:
            c.execute(getSbyUID[self.batch.stage],(self.job_uid,))
            rows=c.fetchall()
            # print(*rows,sep='\n')
        if c.rowcount <1:
            return
        if c.rowcount >1:
            print("Warning, multiple job match, using just the first")
        if self.job_uid is None and rows[0].job_uid is not None:
            print('Warning, failed to pull before executing pull_stats()')
        if self.job_status is None and self.ghost_ref is None:
            self.job_status=1
    def push(self):
        c = None # self.batch.pldb.named_cursor()
        if self.job_uid is None:
            if self.job_status is None:
                self.job_status=0
            if self.batch.pldb.dry_run:
                print('N',end='',file=self.batch.pldb.ww)
            else:
                c = self.batch.pldb.named_cursor()
                c.execute(insJobRef,(self.ghost_ref,self.jobkey.name,self.batch.stage,self.batch.branch,self.batch.name,self.batch.group,self.job_status))
                self.job_uid=c.lastrowid
                self.batch.job_ref[self.jobkey]=self.job_uid
                print('+',end='',file=self.batch.pldb.ww)
        #Consider code to update job status
        
        #Find and update any dangling stats rows
        if self.batch.stage not in updateSJUID.keys():
            return
        if self.batch.pldb.dry_run:
            print('U',end='',file=self.batch.pldb.ww)
        else:
            c = self.batch.pldb.named_cursor()
            c.execute(updateSJUID[self.batch.stage],{'UID':self.job_uid, 'name':self.jobkey.name,'Batch':self.jobkey.Batch})
            if c.rowcount >0 :
                print('^'*c.rowcount,end='',file=self.batch.pldb.ww)
            else:
                print('.',end='',file=self.batch.pldb.ww)
    def proc(self):
        self.pull()
        self.pull_stats()
        self.push()
    def tick(self,text='*'):
        PipelineObj.tick(self,text)
    def push_ticks(self):
        self.pldb.push_ticks()

def get_non_empty_slots(d):
    return { k : v for k,v in d.items() if v is not None }

def items_for_empty_slots(cur,new):
    clr = { k  for k,v in cur.items() if v is None }
    return { k : v for k,v in new.items() if v is not None and k in clr }

# path Row(file_uid=11547, job_uid=5418, file_name=u'Lepthr_coeru_L1985_contigs.fasta', purpose=0, path_uid=24314)
# file Row(file_uid=11547, file_category=u'uce_contigs', file_hash=u'9def5372b13f7e15a8f6ea3327d69604', file_status=1, file_size=1574426, job_uid_temp=4727, file_name_temp=u'Lepthr_coeru_L1985.fasta')

# cache_uid=None,
# file_uid=11547, path_uid=24314, job_uid=5418, file_purpose=0,
# inode=144115271844585995,
# file_type=uce_contigs, file_name=Lepthr_coeru_L1985_contigs.fasta, 
# Batch=stage3-manakin, system=cypress, volume=/lustre, job_name=Lepthr_coeru_L1985
# abspath=/lustre/project/ederrybe/tyranni-uce/stage3-manakin/Lepthr_coeru_L1985/Lepthr_coeru_L1985_contigs.fasta, 

##kd=itertools.izip(*itertools.tee(testF.d().iterkeys()))
##>>> zip(*kd)
# >>> kd=itertools.izip(*itertools.tee(testF.d().iterkeys()))
# >>> tcol,dkey=itertools.izip(*kd)
# >>> ", ".join(tcol)
# 'file_uid, path_uid, file_type, file_name, abspath, Batch, job_uid, system, volume, cache_uid, file_purpose, inode, job_name'
# >>> ", ".join(map(lambda k: "%({})s".format(k),dkey)
# ... )
# '%(file_uid)s, %(path_uid)s, %(file_type)s, %(file_name)s, %(abspath)s, %(Batch)s, %(job_uid)s, %(system)s, %(volume)s, %(cache_uid)s, %(file_purpose)s, %(inode)s, %(job_name)s'
# >>> 



##

def mk_ins(tname,dmap):
    tcol,dkey=itertools.izip(*dmap)
    f=lambda k: "%({})s".format(k)
    return('''INSERT INTO {} ({}) VALUES ({})'''.format(tname,", ".join(tcol),", ".join(map(f,dkey))) )

def mk_upd(tname,dmap,kmap):
    #tcol,dkey=itertools.izip(*dmap)
    f=lambda k: "{}=%({})s".format(k[0],k[1])
    return('''UPDATE {} SET {} WHERE {}'''.format(tname,", ".join(imap(f,dmap))," and ".join(imap(f,kmap)) ))



class JobFile(PipelineObj):
    "Simple class to wrap around a known job file"
    def __init__(self, batch, filekey):
        "Load up parts"
        PipelineObj.__init__(self,batch)
        self.batch_ref   = batch   #Batch object
        self.jobkey  = db2jk(filekey)  #JobKey named tuple
        self.filekey = filekey #FilepathRow named tuple
        self.system=batch.pldb.location
        self.file_type = None
        self.file_purpose = None
        self.file_uid = None
        self.job_uid = None
        self.path_uid = None
        if self.jobkey in batch.job_ref.keys():
            self.job_uid=batch.job_ref[self.jobkey]
        self.cache_uid = None
        self.__dict__.update(filekey.__dict__)
    def d(self):
        "Get dictionary suitable for cursor execution"
        return { k : v for k,v in self.__dict__.items() if k not in ['filekey','jobkey','batch_ref','pldb','auto_commit'] }
    file_uid_row='''SELECT * from files where file_uid=%(file_uid)s'''
    job_uid_row='''SELECT * from jobs where job_uid=%(job_uid)s'''
    job_ref_row='''SELECT * from jobs where job_name=%(job_name)s and batch=%(Batch)s '''
    path_uid_row='''SELECT * from job_files where path_uid=%(path_uid)s'''
    path_ref_row='''SELECT * from job_files where job_uid=%(job_uid)s and file_name=%(file_name)s'''
    path_ref2_row='''SELECT * from job_files where job_uid=%(job_uid)s and file_uid=%(file_uid)s'''
    file_ref_cache='''SELECT distinct file_uid,file_type from file_cache where system=%(system)s and volume=%(volume)s 
    and inode=%(inode)s and file_uid is not null'''
    path_ref_cache='''SELECT * from file_cache where system=%(system)s and abspath=%(abspath)s'''
    path_ref2_cache='''SELECT * from file_cache where system=%(system)s and abspath=%(abspath)s and job_uid=%(job_uid)s'''
    raw_data_row='''SELECT * from raw_data where concat('raw-reads/',rgID,'-READ',R_NUM,'.fastq.gz')=%(file_name)s'''
    cache_data_cols=(u'cache_uid', u'file_uid', u'path_uid', u'job_uid', u'file_purpose', 
    u'file_type', u'inode', u'Batch', u'system', u'volume', u'job_name', u'abspath')
    def data_for_cache(self):
        return { k : v for k,v in self.d().items() if v is not None and k in self.cache_data_cols }
    def load_from_cache(self):
        "Load any info about this file from the cache"
        c = self.batch_ref.pldb.named_cursor()
        # check cache for (system, path) key
        if self.job_uid is not None:
            c.execute(self.path_ref2_cache,self.d())
        else:
            c.execute(self.path_ref_cache,self.d())
        rows = c.fetchall()
        if self.batch_ref.pldb.dry_run: print(*rows,sep='\n')
        if c.rowcount > 1:
            print("Warning, likely issue with cache integrity")
        if c.rowcount >0:
            self.__dict__.update(items_for_empty_slots(self.__dict__,rows[0].__dict__))
            print('V',end='',file=self.batch_ref.pldb.ww)
        # check cache for file_uid and type
        if self.file_uid is None and self.system is not None and self.abspath is not None:
            c.execute(self.file_ref_cache,self.d())
            rows = c.fetchall()
            if c.rowcount > 1:
                print("Warning, likely issue with cache integrity")
            if c.rowcount >0:
                self.__dict__.update(items_for_empty_slots(self.__dict__,rows[0].__dict__))
                print('v',end='',file=self.batch_ref.pldb.ww)
    def check_known_job_files(self):
        "Is this file for this job registered?"
        search = self.path_ref_row
        if self.path_uid is not None:
            if self.file_uid is not None and self.file_purpose is not None:
                return
            search = self.path_uid_row
        elif self.file_uid is not None:
            search = self.path_ref2_row
        c = self.batch_ref.pldb.named_cursor()
        c.execute(search, self.d())
        rows = c.fetchall()
        if self.batch_ref.pldb.dry_run: print(*rows,sep='\n')
        if c.rowcount > 1:
            print("Warning, likely issue with job_files integrity")
        if c.rowcount >0:
            self.file_uid=rows[0].file_uid
            self.path_uid=rows[0].path_uid
            self.file_purpose=rows[0].purpose
            print('p',end='',file=self.batch_ref.pldb.ww)
    def check_known_files(self):
        "Is this file registered?"
        c = self.batch_ref.pldb.named_cursor()
        c.execute(self.file_uid_row,self.d())
        rows = c.fetchall()
        if self.batch_ref.pldb.dry_run: print(*rows,sep='\n')
        if c.rowcount > 1:
            print("Warning, likely issue with files integrity")
        if c.rowcount >0:
            self.file_type=rows[0].file_category
            print('f',end='',file=self.batch_ref.pldb.ww)
    def check_raw_data(self):
        "Look up file in raw data"
        c = self.batch_ref.pldb.named_cursor()
        c.execute(self.raw_data_row,self.d())
        rows = c.fetchall()
        if self.batch_ref.pldb.dry_run: print(*rows,sep='\n')
        if c.rowcount > 1:
            print("Warning, likely issue with raw_data integrity")
        if c.rowcount >0:
            self.file_uid=rows[0].file_uid
            print('r',end='',file=self.batch_ref.pldb.ww)
    def pull(self):
        "Load all available information from database"
        self.load_from_cache()
        self.check_known_job_files()
        if self.file_uid is not None and self.file_type is None:
            self.check_known_files()
        if self.file_purpose is None:
            self.file_purpose=self.det_file_purpose()
        if self.file_type is None:
            self.file_type=self.det_file_type()
        if re.match('raw_read',self.file_type) is not None and self.file_uid is None:
            self.check_raw_data()
        if self.batch_ref.pldb.dry_run: print(self.d())
    ins_file=[('file_category','file_type'),('job_uid_temp','job_uid'),('file_name_temp','file_name')]
    ins_path=[('purpose','file_purpose'),('job_uid','job_uid'),('file_uid','file_uid'),('file_name','file_name')]
    def push_file_row(self,c,all_set):
        "Insert row if needed in files table"
        if self.file_uid is None:
            if self.file_purpose not in [1,3,4] or self.job_uid is None:
                print("CAUTION!",self.file_type,self.file_purpose,self.file_name)
            else:
                if self.batch_ref.pldb.dry_run:
                    print(mk_ins('files',self.ins_file)%all_set)
                else:
                    c.execute(mk_ins('files',self.ins_file),all_set)
                    self.file_uid=c.lastrowid
                    all_set['file_uid']=self.file_uid
                    self.tick('+')
        return(all_set)
    def	push_path_row(self,c,all_set):
       	"Insert	row if needed in job_files table"
        if self.path_uid is None:
            if self.job_uid is None or self.file_uid is None:
                print("UNREGISTERED!",self.job_uid,self.file_uid,self.file_purpose,self.file_name)
            else:
                if self.batch_ref.pldb.dry_run:
                    print(mk_ins('job_files',self.ins_path)%all_set)
                else:
                    c.execute(mk_ins('job_files',self.ins_path),all_set)
                    self.path_uid=c.lastrowid
                    all_set['path_uid']=self.path_uid
                    self.tick('+')
       	return(all_set)
    def push(self):
        "Push collected information to database"
        c = self.batch_ref.pldb.named_cursor()
        all_set=get_non_empty_slots(self.d())
        if self.batch_ref.pldb.dry_run: print(', '.join([ '{0}=%({0})s'.format(k) for k in all_set.keys()]) % all_set)
        if self.job_uid is None: print("ERROR!",self.Batch,self.job_name) 
        if self.file_uid is None: 
            if self.file_purpose not in [1,3,4] or self.job_uid is None:
                print("CAUTION!",self.file_type,self.file_purpose,self.file_name)
            else:
                if self.batch_ref.pldb.dry_run:
                    print(mk_ins('files',self.ins_file)%all_set)
                else:
                    c.execute(mk_ins('files',self.ins_file),all_set)
                    self.file_uid=c.lastrowid
                    all_set['file_uid']=self.file_uid
                    self.tick('+')
        if self.path_uid is None:
            if self.job_uid is None or self.file_uid is None:
                print("UNREGISTERED!",self.job_uid,self.file_uid,self.file_purpose,self.file_name)
            else:
                if self.batch_ref.pldb.dry_run:
                    print(mk_ins('job_files',self.ins_path)%all_set)
                else:
                    c.execute(mk_ins('job_files',self.ins_path),all_set)
                    self.path_uid=c.lastrowid
                    all_set['path_uid']=self.path_uid
                    self.tick('+')
        cache_slots=self.data_for_cache() 
        kd=itertools.izip(*itertools.tee(cache_slots.iterkeys()))
        if self.cache_uid is None:
            if self.batch_ref.pldb.dry_run:
                print(mk_ins('file_cache',kd) % cache_slots)
            else:
                c.execute(mk_ins('file_cache',kd), cache_slots)
                self.cache_uid=c.lastrowid
                self.tick('+')
        else:
            if self.batch_ref.pldb.dry_run:
                print(mk_upd('file_cache',kd,[('cache_uid','cache_uid')]) % cache_slots)
            else:
                c.execute(mk_upd('file_cache',kd,[('cache_uid','cache_uid')]), cache_slots)
                if c.rowcount >0 :
                    self.tick('^'*c.rowcount)
                else:
                    self.tick('.')
    def peek(self):
        cache_row=self.d()
        print(', '.join([ '{0}=%({0})s'.format(k) for k in cache_row.keys()]) % cache_row)
    def det_checksum(self):
        "Run md5sum and return value"
        pass
    def det_file_purpose(self):
        if is_stage1(self.batch_ref.stage):
            if re.match('^raw-data|^raw-reads',self.file_name):
                return 0
            if re.match('^split-adapter',self.file_name):
                return 1
            raise Exception("Unknown file_purpose",self.file_name)
        if is_stage2(self.batch_ref.stage):
            if re.match('^rg....._.._T-READ',self.file_name):
                return 2
            if re.match('^contigs.fa',self.file_name):
                return 4
            if re.search('fasta$',self.file_name):
                return 1
            raise Exception("Unknown file_purpose",self.file_name)
        if is_stage3(self.batch_ref.stage):
            if re.match('^rg....._.._T-READ',self.file_name):
                return 0
            if re.search('_contigs.fasta$',self.file_name):
                return 0
            if re.search('_loci.fasta$',self.file_name):
                return 1
            if re.search('cns.fq.gz$',self.file_name):
                return 3
            if re.search('.depth.txt.gz$',self.file_name):
                return 3
            if re.search('.depth.txt$',self.file_name):
                return 3
            if re.search('.full_depth.gz$',self.file_name):
                return 3
            raise Exception("Unknown file_purpose",self.file_name)
        if is_stage4(self.batch_ref.stage):
            if re.search('.fasta$',self.file_name):
                return 1
            raise Exception("Unknown file_purpose",self.file_name)
        if is_stage5(self.batch_ref.stage):
            if re.search('.fasta$',self.file_name):
                return 1
            raise Exception("Unknown file_purpose",self.file_name)
        raise NotImplementedError("Not yet implemented for given stage",self.batch_ref.stage)
    def det_file_type(self):
        if is_stage1(self.batch_ref.stage):
            if re.match('^raw-data|^raw-reads',self.file_name):
                return u'raw_read'
            if re.match('^split-adapter',self.file_name):
                return u'clean_read'
            raise Exception("Unknown file_type",self.file_name)
        if is_stage2(self.batch_ref.stage):
            if re.match('^rg....._.._T-READ',self.file_name):
                return u'clean_read'
            if re.match('^contigs.fa',self.file_name):
                return u'contigs'
            if re.search('fasta$',self.file_name):
                return u'uce_contigs'
            raise Exception("Unknown file_type",self.file_name)
        if is_stage3(self.batch_ref.stage):
            if re.match('^rg....._.._T-READ',self.file_name):
                return u'clean_read'
            if re.search('_contigs.fasta$',self.file_name):
                return u'uce_contigs'
            if re.search('_loci.fasta$',self.file_name):
                return u'loci_for_sample'
            if re.search('cns.fq.gz$',self.file_name):
                return u'site_quality'
            if re.search('.depth.txt.gz$',self.file_name):
                return u'depth_distribution'
            if re.search('.depth.txt$',self.file_name):
                return u'depth_distribution'
            if re.search('.full_depth.gz$',self.file_name):
                return u'site_depth'
            raise Exception("Unknown file_type",self.file_name)
        if is_stage4(self.batch_ref.stage):
            if re.search('fasta$',self.file_name):
                return u'samples_for_locus'
            raise Exception("Unknown file_type",self.file_name)
        if is_stage5(self.batch_ref.stage):
            if re.search('fasta$',self.file_name):
                return u'aligned_sequences'
            raise Exception("Unknown file_type",self.file_name)
        raise NotImplementedError("Not yet implemented for given stage",self.batch_ref.stage)

class BatchJobFiles(PipelineObj):
    "Use a temporary table to register multiple JobFiles"
    def __init__(self, batch, drop_prior = True):
        "Load up parts"
        PipelineObj.__init__(self,batch)
        self.batch_ref = batch   #Batch object
        # self.table_key = table_key

    key_ref_cols=(u'temp_uid', u'file_uid', u'path_uid', u'job_uid', u'inode', u'Batch', u'system', u'volume', u'job_name', u'abspath')
    drop_temp = '''DROP TEMPORARY TABLE IF EXISTS bjf_temp'''
    make_temp = '''CREATE TEMPORARY TABLE IF NOT EXISTS bjf_temp (temp_uid int(11) PRIMARY KEY,
        cache_uid int(11), file_uid int(11), path_uid int(11), job_uid int(11), inode bigint(20),
        Batch varchar(45), system varchar(45), volume varchar(45), file_name varchar(255),
        job_name varchar(45), file_type varchar(45), file_purpose varchar(45), abspath varchar(511))'''
    def d(self):
        "Get dictionary suitable for cursor execution"
        return { k : v for k,v in self.__dict__.items() if k not in ['filekey','jobkey','batch_ref','table_ref', 'pldb','auto_commit'] }
    file_purp_s45='''UPDATE bjf_temp t set file_purpose=1 where file_name regexp 'fasta$' and file_purpose is null and stage regexp '^stage[45]' '''
    file_purp_s3_0='''UPDATE bjf_temp t set file_purpose=0 where file_type in ('clean_reads', 'uce_contigs') and file_purpose is null and stage = 'stage3' '''
    file_purp_s2_0='''UPDATE bjf_temp t set file_purpose=2 where file_type = 'clean_reads' and file_purpose is null and stage = 'stage2' '''
    file_purp_s1_1='''UPDATE bjf_temp t set file_purpose=1 where file_type = 'clean_reads' and file_purpose is null and stage = 'stage1' '''
    file_purp_s1_0='''UPDATE bjf_temp t set file_purpose=0 where file_type = 'raw_reads' and file_purpose is null and stage = 'stage1' '''
    file_type_s5='''UPDATE bjf_temp t set file_type='aligned_sequences' where file_name regexp 'fasta$' and file_type is null'''
    file_type_s4='''UPDATE bjf_temp t set file_type='samples_for_locus' where file_name regexp 'fasta$' and file_type is null'''
    file_type_s1_A='''UPDATE bjf_temp t set file_type='raw_read' where file_name regexp '^raw' and file_type is null'''
    file_type_s1_B='''UPDATE bjf_temp t set file_type='clean_read' where file_name regexp '^split-adapter' and file_type is null'''
    file_type_s2_A='''UPDATE bjf_temp t set file_type='clean_read' where file_name regexp '^rg....._.._T-READ' and file_type is null'''
    file_type_s2_B='''UPDATE bjf_temp t set file_type='contigs' where file_name regexp '^contigs.fa' and file_type is null'''
    file_type_s2_C='''UPDATE bjf_temp t set file_type='uce_contigs' where file_name regexp 'fasta$' and file_type is null'''
    file_type_s3_A='''UPDATE bjf_temp t set file_type='clean_read' where file_name regexp '^rg....._.._T-READ' and file_type is null'''
    file_type_s3_B='''UPDATE bjf_temp t set file_type='uce_contigs' where file_name regexp 'contigs.fasta$' and file_type is null'''
    file_type_s3_C='''UPDATE bjf_temp t set file_type='loci_for_sample' where file_name regexp 'loci.fasta$' and file_type is null'''
    file_type_s3_D='''UPDATE bjf_temp t set file_type='site_quality' where file_name regexp 'cns.fq.gz$' and file_type is null'''
    file_type_s3_E='''UPDATE bjf_temp t set file_type='depth_distribution' where file_name regexp '.depth.txt.gz$' and file_type is null'''
    file_type_s3_F='''UPDATE bjf_temp t set file_type='depth_distribution' where file_name regexp '.depth.txt$' and file_type is null'''
    file_type_s3_G='''UPDATE bjf_temp t set file_type='site_depth' where file_name regexp '.full_depth.gz$' and file_type is null'''
    file_uid_row='''UPDATE files f join bjf_temp t using(file_uid)
        set t.file_type=f.file_category where t.file_type is null'''
    path_uid_row='''update bjf_temp t join job_files f using(path_uid)
        set t.file_uid=f.file_uid, t.file_purpose=f.purpose
        where t.file_uid is null or t.file_purpose is null'''
    path_ref_row='''update bjf_temp t join job_files f using(job_uid,file_name)
        set t.path_uid=f.path_uid, t.file_purpose=f.purpose, t.file_uid=f.file_uid 
        where t.path_uid is null and (t.file_purpose is null or f.file_uid is null)'''
    path_ref2_row='''update bjf_temp t join job_files f using(job_uid,file_uid)
        set t.path_uid=f.path_uid, t.file_purpose=f.purpose 
        where t.path_uid is null and t.file_purpose is null'''
    file_ref_cache='''UPDATE file_cache f join bjf_temp t using(system,volume,inode) 
        SET t.file_uid=f.file_uid,t.file_purpose=f.file_purpose
        where f.file_uid is not null and t.cache_uid is null and (t.file_uid is null or t.file_type is null)'''
    path_ref_cache='''UPDATE file_cache f join bjf_temp t using(system,abspath) 
        SET t.file_uid=ifnull(t.file_uid,f.file_uid),
        t.file_purpose=ifnull(t.file_purpose,f.file_purpose),
        t.job_uid=ifnull(t.job_uid,f.job_uid),
        t.path_uid=ifnull(t.path_uid,f.path_uid),
        t.cache_uid=ifnull(t.cache_uid,f.cache_uid),
        t.file_type=ifnull(t.file_type,f.file_type)
        where t.job_uid is null'''
    path_ref2_cache='''UPDATE file_cache f join bjf_temp t using(system,abspath,job_uid) 
        SET t.file_uid=ifnull(t.file_uid,f.file_uid),
        t.file_purpose=ifnull(t.file_purpose,f.file_purpose),
        t.path_uid=ifnull(t.path_uid,f.path_uid),
        t.cache_uid=ifnull(t.cache_uid,f.cache_uid),
        t.file_type=ifnull(t.file_type,f.file_type)
        '''
    raw_data_row='''UPDATE raw_data f join bjf_temp t 
        on concat('raw-reads/',f.rgID,'-READ',f.R_NUM,'.fastq.gz')=t.file_name
        SET t.file_uid=f.file_uid
        where t.file_uid is null and t.file_type = 'raw-read' '''
    cache_data_cols=(u'cache_uid', u'file_uid', u'path_uid', u'job_uid', u'file_purpose',
    u'file_type', u'inode', u'Batch', u'system', u'volume', u'job_name', u'abspath')

def process_batch_structure(bT):
    for j in bT.iter_batch_jobs(): j.proc();
    bT.flush_ticks()
    tmp=[(jF.pull(),jF.push()) for jF in bT.iter_job_files()]
    print(len(tmp))
    bT.push_ticks()



def process_batch_specials(bT, job_status):
    for j in bT.iter_batch_jobs(): 
        j.pull()
        if j.job_uid is None:
            j.job_status=job_status
        j.push()
    bT.push_ticks()


job_header=["Branch","BatchGroup","Batch","SampleName"]
s3_stat=["nLoci","avg_length","avg_coverage","avg_depth"]
stat_footer=["Order_no","sample_id","job_uid","file_exists","file_accepted"]
Stage3Stat=namedtuple("Stage3Stat",job_header+["LociFileName"]+s3_stat+stat_footer)
empty_s3_stat=Stage3Stat(**dict(zip(Stage3Stat._fields,itertools.repeat(None))))
ins_s3_stat=mk_ins('stage3stats',zip(Stage3Stat._fields,Stage3Stat._fields))
upd_s3_stat=mk_upd('stage3stats',zip(Stage3Stat._fields,Stage3Stat._fields),[('job_uid','job_uid')])

def default_stage3_stats(job):
    return empty_s3_stat._replace(job_uid=job.job_uid,Branch=job.batch.branch,
        Batch=job.batch.name,BatchGroup=job.batch.group,SampleName=job.jobkey.name,
        LociFileName=job.jobkey.name+"_loci.fasta" )

## test code
# python -i load_tyranniSQL.py 
# bT=Batch(pdb,batchpath="stage3-Thamno_nigro");jFs=bT.iter_job_files()
# jF=next(jFs); jF.pull()
# depth_file=tyranniSQL.PipelineDb.FullDepthFile(jF)
# depth_file.pull()

# python -i load_tyranniSQL.py 
# bT=Batch(pdb,batchpath="stage3-Thamno_nigro");
# for jF in bT.iter_job_files():
#     jF.pull()
#     jF.push()
#     if jF.file_type == u'site_depth':
#         depth_file=tyranniSQL.PipelineDb.FullDepthFile(jF)
#         depth_file.pull()
#         depth_file.parse_file()


class FullDepthFile(JobFile):
    "Ecapsulates a full depth file"
    def __init__(self, batch, filekey=None):
        if issubclass(type(batch),JobFile):
            #Upgrade job file
            JobFile.__init__(self,batch=batch.batch_ref,filekey=batch.filekey)
        else:
            JobFile.__init__(self,batch=batch,filekey=filekey)
        self.ins_depth_stmt =  """insert into sparse_sequence_depth 
        (job_uid,locus_uid,seq_uid, site, depth) VALUES
        ({},?,?,?,?)""".format(self.job_uid)
    def pull(self,upload_all=False):
        JobFile.pull(self)
        if self.file_type != u'site_depth':
            raise Exception("File is not a site depth file!")
        if self.file_uid is None: 
            print("Warning, unregistered file!")
            JobFile.push(self)
        c=self.named_cursor()
        c.execute("""SELECT file_uid fasta_file
        FROM  job_files where job_uid=%s and purpose=1""",(self.job_uid,))
        self.fasta_file_uid=c.fetchall()[0].fasta_file
#         if upload_all:
#             c.execute("""SELECT locus_name, seq_uid 
#             FROM  job_files natural join unaligned_sequences 
#                 natural join loci where job_uid=%s""",(self.job_uid,))
#         else:
#             c.execute("""SELECT locus_name, seq_uid 
#             FROM  job_files natural join unaligned_sequences 
#                 natural join loci where job_uid=%s and depth_data is null""",(self.job_uid,))
        if upload_all:
            c.execute("""SELECT locus_uid, seq_uid, length(seq_data) len 
            FROM unaligned_sequences where file_uid=%s 
            for update""",(self.fasta_file_uid,))
        else:
            c.execute("""SELECT locus_uid, seq_uid, length(seq_data) len 
            FROM unaligned_sequences where file_uid=%s and depth_data is null
            for update""",(self.fasta_file_uid,))
        rows=[(str(self.pldb.get_locus(k)),v,l) for k,v,l in c]
        self.seq_d=dict(((k,v) for k,v,l in rows))
        self.seqs=[ v for k,v,l in rows if v is not None]
        self.loci=[ k for k,v,l in rows if v is not None]
        self.seq_len=dict(((k,l) for k,v,l in rows))
#         self.seq_len=dict()
#         for r in self.seq_d.items():
#         c.executemany("""SELECT ? locus_name, length(seq_data) seq_len 
#             FROM  unaligned_sequences where seq_uid=?""",self.seq_d.items())
#         c.execute("""SELECT locus_name, length(seq_data) len 
#         FROM  job_files natural join unaligned_sequences 
#             natural join loci where job_uid=%s""",(self.job_uid,))
#         self.seq_len=dict(((str(k),v) for k,v in c.fetchall()))
    def get_seq(self,locus):
        if locus not in self.loci:
            return SequenceDepth(p_obj=self,seq_uid=None,locus_name=locus)
        return SequenceDepth(p_obj=self,seq_uid=self.seq_d[locus],locus_name=locus,n_sites=self.seq_len[locus])
    def parse_file(self,post_large_depths=False,only_large_depths=False):
        with io.TextIOWrapper(io.BufferedReader(gzip.open(self.abspath))) as src:
            lines = iter(src)
            l = next(lines)
            v=reD.match(l).groups()
            cur_seq=self.get_seq(v[0])
            cur_seq[int(v[1])]=int(v[2])
            if self.pldb.dry_run:
                print(v[0],cur_seq._name,end='')
                self.pldb.ww.write('>')
            for l in lines:
                #self.tick('-')
                v=reD.match(l).groups()
                if cur_seq._name != v[0]:
                    if self.pldb.dry_run:
                        self.pldb.ww.write('/')
                    cur_seq.push(post_large_depths=post_large_depths,only_large_depths=only_large_depths)
                    if not self.pldb.rollback:
                        self.pldb.conn.commit()
                    cur_seq=self.get_seq(v[0])
                #print(v,end='') #find that bug!
                cur_seq[int(v[1])]=int(v[2])
            if self.pldb.dry_run:
                self.pldb.ww.write(':')
        cur_seq.push(post_large_depths=post_large_depths,only_large_depths=only_large_depths)
        if not self.pldb.rollback:
            self.pldb.conn.commit()
        self.pldb.push_ticks()
    def parse_file_burst(self,post_large_depths=False,only_large_depths=False):
        with io.TextIOWrapper(io.BufferedReader(gzip.open(self.abspath))) as src:
            lines = iter(src)
            l = next(lines)
            v=reD.match(l).groups()
            read_seq=list()
            cur_seq=self.get_seq(v[0])
            cur_seq[int(v[1])]=int(v[2])
            if self.pldb.dry_run:
                print(v[0],cur_seq._name,end='')
                self.pldb.ww.write('>')
            for l in lines:
                #self.tick('-')
                v=reD.match(l).groups()
                if cur_seq._name != v[0]:
                    if self.pldb.dry_run:
                        self.pldb.ww.write('/')
                    if not only_large_depths:
                        cur_seq.prep()
                    read_seq.append(cur_seq)
                    cur_seq=self.get_seq(v[0])
                #print(v,end='') #find that bug!
                cur_seq[int(v[1])]=int(v[2])
        if self.pldb.dry_run:
            self.pldb.ww.write(':')
        for s in read_seq:
            s.push(post_large_depths=post_large_depths,only_large_depths=only_large_depths)
            if not self.pldb.rollback:
                self.pldb.conn.commit()
        cur_seq.push(post_large_depths=post_large_depths,only_large_depths=only_large_depths)
        if not self.pldb.rollback:
            self.pldb.conn.commit()
        self.pldb.push_ticks()

SiteDepth = namedtuple("SiteDepth",["locus","site","depth"])
FilePathTag = namedtuple("FilePathTag",["batch","job_name","file_name"])
FileUIDTag = namedtuple("FileUIDTag",["job_uid","file_uid","path_uid","file_path"])


sql_sparse_depth = """
create temporary table sparse_sequence_depth
(seq_uid int,locus_uid int, site int, depth int);"""

#>>> lines = gzip.GzipFile("stage3-282462/Myiorn_aur_MZU79945/Myiorn_aur_MZU79945.full_depth.gz","r")
#>>> l = next(lines)
reD = re.compile("^([^|]+)\|[^|]+\t([0-9]+)\t([0-9]+)")





class SequenceDepth(PipelineObj):
    def __init__(self,p_obj,seq_uid,locus_name,n_sites=None):
        PipelineObj.__init__(self,p_obj)
        self.seq_uid=seq_uid
        self._name=locus_name
        self.locus_uid=self.pldb.get_locus_uid(locus_name)
        if self.seq_uid is not None:
            if n_sites is None:
                c=self.named_cursor()
                c.execute("""SELECT length(seq_data) len
                FROM unaligned_sequences
                where seq_uid=%s""",(self.seq_uid,))
                n_sites=c.fetchall()[0][0]
            self._data=bytearray(n_sites)
        else:
            self._data=None
            self.pldb.ww.write('*')
        self._sparse=dict()
        self.sub_start=None
        self.sub_end=None
        self.sub_avg=None
        self.seq_avg=None
        self.seq_qual=None
        if self.pldb.dry_run:
            self.pldb.ww.write('|')
    def push(self,post_large_depths=False,only_large_depths=False):
        if self.seq_uid is None: return
        c=self.named_cursor()
        compact_args=None
        if not only_large_depths:
            compact_args=dict()
            compact_args['seq_uid']=self.seq_uid
            compact_args['data']=self._data
            if self.seq_qual is None:
                compact_args['hi_qual']=int(self.check_quality())
            else: compact_args['hi_qual']=int(self.seq_qual)
            compact_args['sub_start']=self.sub_start
            compact_args['sub_end']=self.sub_end
            if self.seq_avg is None:
                compact_args['avg_depth']=self.get_avg()
            else: compact_args['avg_depth']=self.seq_avg
            if self.sub_avg is None:
                compact_args['sub_depth']=self.get_sub_avg()
            else: compact_args['sub_depth']=self.sub_avg
        if post_large_depths:
            sp_data=[ (self.seq_uid,self.locus_uid,site,depth) for site,depth in self._sparse.iteritems()]
        if self.pldb.dry_run:
            #print(self._name)
            if compact_args is None: self.tick('?')
            else: self.tick('U')
            #for i,d in itertools.izip(range(3),sp_data):
            #    self.tick('^')
            return
        marks=0;
        if not only_large_depths:
            c.execute("""UPDATE unaligned_sequences set
            depth_data=%(data)s,
            sub_start =%(sub_start)s,
            sub_end   =%(sub_end)s,
            avg_depth =%(avg_depth)s,
            sub_depth =%(sub_depth)s,
            hi_qual =%(hi_qual)s
            where seq_uid=%(seq_uid)s""",compact_args)
            marks+=c.rowcount
        ## iterate singles...?"
        if  post_large_depths:
            if len(sp_data) > 0:
                c.executemany("""insert into sparse_sequence_depth (seq_uid,locus_uid,site,depth)
        VALUES (%s,%s,%s,%s) ON DUPLICATE KEY UPDATE depth = VALUES(depth)""",sp_data)
                marks+=c.rowcount
        if marks>0: self.tick('+')
        else: self.tick('-')
    def prep(self):
        self.seq_avg=self.get_avg()
        self.seq_qual=self.check_quality()
        self.sub_avg=self.get_sub_avg()
    def __getitem__(self,i):
        v=self._data[i-1]
        if v is 255:
            return self._sparse[i]
        else:
            return v
    def __setitem__(self,i,value):
        if self.seq_uid is None: return
        if i < 1:
            self.tick('!<')
            if self.pldb.dry_run: print(self._name,i,value,end='')
            return
        if i > len(self._data):
            if value < 10: return ##Ignore bad index for low depth sites
            self.tick('!>') ##Issue an alert, but keep going; it will not cause issues
            if self.pldb.dry_run: print(self._name,i,value,end='')
            return
        if value >= 255:
            self._data[i-1]=255
            self._sparse[i]=value
        else:
            self._data[i-1]=value
            if i in self._sparse.keys():
                self._sparse.pop(i)
        if value >=5:
            self.sub_end=i
            if self.sub_start is None:
                self.sub_start = i-1
    def get_sub_avg(self):
        if self.sub_start is None: return 0;
        if self.sub_end is None: return 0;
        res=sum(self._data[self.sub_start:self.sub_end])
        for k,v in self._sparse.iteritems():
            if k <= self.sub_start or k > self.sub_end:
                continue
            res+=v-255
        return float(res)/float(self.sub_end-self.sub_start)
    def get_avg(self):
        if self.seq_uid is None: 
            self.tick('*')
            return 0
        res=sum(self._data)
        for k,v in self._sparse.iteritems():
            res+=v-255
        return float(res)/float(len(self._data))
    def check_quality(self):
        #qual=False
        if self.sub_start is None:
            return False
        if self.sub_end-self.sub_start <=100:
            return False
        if self.get_sub_avg() <15:
            while(self.shrink_sequence()):
                if self.get_sub_avg() >=15:
                    return True
        else: return True
        return False
    def shrink_sequence(self):
        if self.sub_start is None: return False
        if self.sub_end-self.sub_start <=100:
            return False
        if self._data[self.sub_end-1] > b'\x0E' and self._data[self.sub_start] > b'\x0E':
            return False
        if self._data[self.sub_end-1] == self._data[self.sub_start]:
            self.sub_end-=1
            if self.sub_end-self.sub_start >100:
                self.sub_start+=1
            return True
        if self._data[self.sub_end-1] < self._data[self.sub_start]:
            self.sub_end-=1
        else:
            self.sub_start+=1
        return True

def collect_stage4_records(jobfile):
    "Upload the record keys from the stage4 file to the database"
    if re.match('^stage4',jobfile.batch_ref.stage) is None:
        raise Exception("Can only collect stage4 records from stage4 files!")
    # Make sure job is registered
    if jobfile.file_uid is None:
        jobfile.pull()
        jobfile.push()
    # Check to see if stats are already present
    c=jobfile.named_cursor()
    c.execute('''SELECT count(*) seqs from stage4sequences where file_uid =%s limit 5''',(jobfile.file_uid,))
    ## is this job already listed? 
    rows=c.fetchall()
    if rows[0].seqs >0:
        jobfile.tick('.'*rows[0].seqs)
        jobfile.pldb.push_ticks()
        return(True)
    filepath=jobfile.abspath
    insrec="""INSERT IGNORE into stage4sequences (seq_uid,sample_uid,locus_uid,file_uid)
    VALUES (%(SEQUID)s,%(SUID)s,%(LUID)s,{})""".format(jobfile.file_uid)
    if os.stat(filepath)[6] > 3:
        for seq_record in SeqIO.parse(filepath, "fasta"):
            parts = seq_record.id.split('|')
            if len(parts) < 5:
                raise Exception('Missing metadata',seq_record.id, filepath)
            rec=dict([tuple(i.split('_')) for i in parts[2:5]])
            if jobfile.pldb.dry_run:
                jobfile.tick('U')
                continue
            c.execute(insrec,rec)
            if c.rowcount >0:
                jobfile.tick('+')
            else:
                jobfile.tick('-')
        jobfile.pldb.push_ticks()
        return(True)
    else:
        jobfile.tick('!')
        return(False)

def collect_stage3_stats(job):
    "Using the record of the job from the database, post the stage3 stats"
    if re.match('^stage3',job.batch.stage) is None:
        raise Exception("Can only collect stage3 stats from stage3 jobs!")
    # Make sure job is registered
    if job.job_uid is None:
        job.proc()
    # Check to see if stats are already present
    c=job.named_cursor()
    c.execute('''SELECT * from stage3stats where job_uid =%s''',(job.job_uid,))
    rows=c.fetchall()
    ## is this job already listed? 
    if c.rowcount >0:
        # Assume only 1 row
        oldrow=rows[0]
        modrow=default_stage3_stats(job)._replace(**get_non_empty_slots(oldrow._asdict()))
        loader=modrow._asdict()
        #print("Warning: not checking to see if anything needs to be fixed!")
        job.tick('.')
        return(True)
    job_path=os.path.join(job.batch.path,job.jobkey.name)
    files=os.listdir(job_path)
    if job.jobkey.name+".depth.txt" not in files and job.jobkey.name+"_loci.fasta" not in files:
        job.tick('-')
        return(False)
    oldrow= default_stage3_stats(job)._replace(file_exists=1,file_accepted=0)
    loader=oldrow._asdict()
    #Look up SampleName and get ID codes
    c.execute('''SELECT sample_uid,alias_uid from sample_aliases where SampleName=%(SampleName)s''',loader)
    id=c.fetchall()[0]
    loader['sample_id']=id.alias_uid
    loader['Order_no']=id.sample_uid
    if oldrow.SampleName+".depth.txt" in files:
        loader['avg_depth']=parse_depth_file(os.path.join(job_path,oldrow.SampleName+".depth.txt"))
    if oldrow.SampleName+"_loci.fasta" in files:
        loader['nLoci'],loader['avg_length'],loader['avg_coverage']=parse_loci_file(os.path.join(job_path,oldrow.SampleName+"_loci.fasta"))
    if job.batch.pldb.dry_run:
        job.tick('N')
        print(Stage3Stat(**loader))
        return(True)
    else:
        c.execute(ins_s3_stat,loader)
        if c.rowcount >0:
            job.tick('+')
            return(True)
        else:
            job.tick('.')
            return(False)
    #raise NotImplementedError("Collecting stage3 stats not yet implemented")

s2sRE=re.compile('contigs:')
s2pRE=re.compile('''([0-9]+) \(([0-9]+\.[0-9]+)%\) uniques of ([0-9]+) contigs,'''+
''' ([0-9]+) dupe probe matches, ([0-9]+) [^,]*, ([0-9]+)''')

s2head='''SampleName,Branch,BatchGroup,Batch,unique_contigs,percent_unique,contigs,
    dupe_probe,multi_contig_matches,multi_probe_matches,
    sample_id,order_no,job_uid,file_exists,file_accepted'''
s2row='''%(name)s,%(branch)s,%(batchgroup)s,%(batch_name)s,
    %(unique)s,%(p_unique)s,%(contigs)s,%(dupes)s,%(multi_probe)s,%(multi_contig)s,
    %(sample_id)s,%(Order_no)s,%(job_uid)s,1,0'''
stage2ins='''INSERT INTO stage2stats ({0}) VALUES ({1})'''.format(s2head,s2row)
s2stat=['unique','p_unique','contigs','dupes','multi_probe','multi_contig']
#s2uid=['sample_id','order_no','job_uid']

def collect_stage2_stats(job):
    "Using the record of the job from the database, post the stage2 stats"
    if re.match('^stage2',job.batch.stage) is None:
        raise Exception("Can only collect stage2 stats from stage2 jobs!")
    # Make sure job is registered
    if job.job_uid is None:
        job.proc()
    # Check to see if stats are already present
    c=job.named_cursor()
    c.execute('''SELECT count(*) n from stage2stats where job_uid =%s''',(job.job_uid,))
    rows=c.fetchall()
    ## is this job already listed? 
    if rows[0].n >0:
        # Assume only 1 row, data collected, skip stats for this job
        # oldrow=rows[0]
        #modrow=default_stage2_stats(job)._replace(**get_non_empty_slots(oldrow._asdict()))
        #loader=modrow._asdict()
        #print("Warning: not checking to see if anything needs to be fixed!")
        job.tick('.')
        return(True)
    job_path=os.path.join(job.batch.path,job.jobkey.name)
    files=os.listdir(job_path)
    if "extract_uce_log.txt" not in files:
        job.tick('-')
        return(False)
    stats={'batch_name':job.batch.name,'branch':job.batch.branch,
        'batchgroup':job.batch.group}
    stats.update(job.jobkey._asdict())
    c.execute('''SELECT sample_uid,alias_uid from sample_aliases 
    where SampleName=%s''',(job.jobkey.name,))
    row=c.fetchall()[0]
    stats['Order_no']=row.sample_uid
    stats['sample_id']=row.alias_uid
    stats['job_uid']=job.job_uid
    with open(os.path.join(job_path,"extract_uce_log.txt")) as f:
        for line in f:
            if s2sRE.search(line) is not None:
                res=s2pRE.search(line)
                if res is not None:
                    stats.update(zip(s2stat,res.groups()))
                    break
    if job.batch.pldb.dry_run:
        job.tick('N')
        print(stats)
        return(True)
    else:
        c.execute(stage2ins,stats)
        if c.rowcount >0:
            job.tick('+')
            return(True)
        else:
            job.tick('.')
            return(False)


def collect_stage1_stats(job):
    "Using the record of the job from the database, post the stage1 stats"
    if re.match('^stage1',job.batch.stage) is None:
        raise Exception("Can only collect stage1 stats from stage1 jobs!")
    raise NotImplementedError("Collecting stage1 stats not yet implemented")


def parse_loci_file(filepath):
    "Calculates number of loci, mean length and mean coverage for the given fasta file"
    if os.stat(filepath)[6] > 3:    
        lengths = list()
        coverages = list()
        i = 0
        for seq_record in SeqIO.parse(filepath, "fasta"):
            parts = seq_record.id.split('|')
            if len(parts) < 2:
                raise Exception('Missing metadata',seq_record.id, filepath)
            nodename = parts[1]
            nodename_parts = nodename.split('_')
            coverages.append(float(nodename_parts[5]))
            lengths.append(len(seq_record))
            i += 1
        if len(lengths) > 0:
            avg_length = float(sum(lengths)/len(lengths))
        else:
            avg_length = 0
        if len(coverages) > 0:
            avg_coverage = float(sum(coverages)/len(coverages))
        else:
            avg_coverage = 0
        return( i, avg_length, avg_coverage)
    else:
        return (None,None,None)


def parse_loci_file_db(file_record):
    "Calculates number of loci, mean length and mean coverage for the given file in the db"
    raise NotImplementedError("Collecting stage3 stats from db not yet implemented")

def parse_depth_file(filepath):
    "Calculates the mean depth from depth log file"
    if os.stat(filepath)[6] > 0:
        infile = open(filepath,"rU")
        depth_counts = list()
        counts = list()
        for line in infile:
            parts = line.split()
            count = int(parts[0])
            counts.append(count)
            depth = int(parts[1])
            depth_counts.append(count*depth)
        return(float(sum(depth_counts))/sum(counts))
    else:
        return(None)
