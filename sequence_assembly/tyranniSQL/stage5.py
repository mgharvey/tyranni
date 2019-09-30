#!/usr/bin/env python2
##useful

from __future__ import print_function
import itertools, os
import datetime,textwrap
from Bio import SeqIO
import mysql.connector
import re,os.path,sys,string
#import openpyxl,lxml
#from openpyxl import load_workbook
from mysql.connector.errors import IntegrityError
"""

Name: stage5.py

Author: Graham E Derryberry
Date: 13 September 2016

Description:
Interface stage 5 fasta files with the mysql db.


"""

class Stage5SamplesFile:
    """Represents a stage 5 fasta file in the database"""
    def __init__(self, pldb,fp=None, locus_uid=None,locus_name=None,group_name=None,delay=False):
        if delay:
            raise NotImplementedError("Delay connection code not yet implemented")
        else:
            self.pldb=pldb
        self.pldb.ww.tick=50
        self.locus_name =locus_name
        self.locus_uid=locus_uid
        self.filepath=fp
        self.group=group_name
        self.aln_len=None
#         if self.locus_uid is not None:
#             self.fetch_rows=''' select  sample_uid,seq_uid,locus_uid,aln_uid,SampleName,
#             aligned_seq_data from aligned_sequences join sample_aliases using(sample_uid)
#             where prefer_align=b'1' and locus_uid={} and file_uid={}'''.format(self.locus_uid,self.file_uid)
#         else:
        self.fetch_rows=None
        if pldb.conn is not None and self.filepath is not None:
            self.load_file_ids()
        elif pldb.conn is not None and self.group is not None:
            self.select_file_uid()
        else:
            self.file_uid=None
        print('V',end='',file=sys.stderr)
    select_src_id='''select 
        file_uid, Order_no as sample_uid, sample_id as alias_uid, branch 
    from (select file_uid, job_uid, eff_job_id from known_output_files 
        where file_category='loci_for_sample' and file_path=%s ) as src
    join stage3stats as j on j.job_uid=src.job_uid or j.job_uid=src.eff_job_id'''
    insRow = '''insert ignore into aligned_sequences 
        (sample_uid,locus_uid,reversed,aligned_seq_data,seq_uid,file_uid) values 
        (?,?,?,?,?,?)'''
    insLen = '''insert ignore into alignment_length
        (file_uid,len) values 
        (?,?)'''
    def select_file_uid(self):
        c=self.pldb.named_cursor()
        if self.locus_uid is not None:
            c.execute('''SELECT aln_file_uid from stage5jobs join jobs using(job_uid) 
                where batchgroup=%s and locus_uid=%s''',(self.group,self.locus_uid))
        elif self.locus_name is not None:
            c.execute('''SELECT aln_file_uid,locus_uid from stage5jobs join jobs using(job_uid) 
                where batchgroup=%s and job_name=%s''',(self.group,self.locus_name))
        else:
            raise NotImplementedError('Sequence download of unknown loci not implemented')
        rows=c.fetchall()
        if c.rowcount is 0:
            raise Exception("Alignment file not found",self.locus_name,self.group)
        self.file_uid=rows[0].aln_file_uid
        if self.locus_uid is None:
            self.locus_uid=rows[0].locus_uid
        self.fetch_rows=''' select  sample_uid,seq_uid,locus_uid,aln_uid,SampleName,
            aligned_seq_data from aligned_sequences join sample_aliases using(sample_uid)
            where prefer_align=b'1' and locus_uid={} and file_uid={}'''.format(self.locus_uid,self.file_uid)
    def load_file_ids(self):
        #self.job_file=self.pldb.load_job_file_from_cache(self.filepath)
        self.file_uid=self.pldb.load_file_uid_from_cache(self.filepath)
#         if self.pldb.conn is None:
#             self.pldb.open_connection()
#         "NEEDED? Get sample uid to sample name for stage 4 alignment dictionary"
#         if self.pldb.conn is None:
#             self.pldb.open_connection()
#         self.file_ids=self.pldb.get_sample_dict(stage='stage4')
#         if self.fetch_rows is None:
#             self.fetch_rows=''' select  sample_uid,branch,seq_data
#             from best_unaligned_sequences join unaligned_sequences using (seq_uid,sample_uid,locus_uid)
#             where locus_uid={}'''.format(self.pldb.get_locus_uid(self.locus_name))
#             self.locus_uid=self.pldb.get_locus_uid(self.locus_name)
    row_header=""">{SampleName}|ALNUID_{aln_uid}|SEQUID_{seq_uid}|SUID_{sample_uid}|LUID_{locus_uid}"""
    def download_sequences(self,out=sys.stdout,filter=None):
        #raise NotImplementedError('Sequence download not yet implented')
        if self.file_uid is None:
            self.select_file_uid()
        c=self.pldb.named_cursor()
        #print(self.insRow);
        #raise NotImplementedError('Sequence upload not yet implented')
        #raise NotImplementedError('Sequence scan not yet implented')
        c.execute(self.fetch_rows)
        with self.pldb.ww.lock:
            if not self.pldb.ww.ticking:
                self.pldb.ww.resume_place()
        tick_cycle = itertools.cycle( ( i % 100 == 0 for i in itertools.count(1) ) )
        if filter is None:
            for row in c:
                print(self.row_header.format(**row._asdict()),file=out)
                print(*wrap(str(row.aligned_seq_data),70),sep="\n",file=out)
                #print(*textwrap.wrap(str(row.aligned_seq_data),70),sep="\n",file=out)
                if next( tick_cycle):
                    print('+',end='',file=self.pldb.ww)
        else:
            for row in c:
                if not filter(row):
                    print('S',end='',file=self.pldb.ww)
                    continue
                print(self.row_header.format(**row._asdict()),file=out)
                print(*wrap(str(row.aligned_seq_data),70),sep="\n",file=out)
                #print(*textwrap.wrap(str(row.aligned_seq_data),70),sep="\n",file=out)
                if next( tick_cycle):
                    print('+',end='',file=self.pldb.ww)
    def build_row(self,SLI,rvd_seq,seq):
        if 'SEQUID' not in SLI.keys():
            raise NotImplementedError("Reverse lookup of seq_uid from src files not implemented")
        return((int(SLI['SUID']),int(SLI['LUID']),rvd_seq,format(seq),int(SLI['SEQUID']),self.file_uid))
    def scan_sequences(self):
        #raise NotImplementedError('Sequence scan not yet implented')
        #if self.file_ids is None:
        if self.pldb.conn is None:
            self.pldb.open_connection()
            #self.load_file_ids()
        print(self.insRow);
        #raise NotImplementedError('Sequence scan not yet implented')
        if os.stat(self.filepath)[6] > 3:
            for i,seq_record in itertools.izip(range(3),SeqIO.parse(self.filepath, "fasta")):
                parts = seq_record.id.split('|')
                if len(parts) < 3:
                    raise Exception('Missing metadata',seq_record.id, self.filepath)
                #node_id = self.pldb.get_locus_uid(parts[0])
                reversed=re.match('^_R_',parts[0]) is not None
                req_id = dict(map(lambda x: tuple(x.split('_')),parts[2:]))
                #s=dict(zip(stats[0::2],stats[1::2]))
                print(self.build_row(req_id,reversed,seq_record.seq))
                # print(dict(zip(stats[0::2],stats[1::2])))
                # print(seq_record.seq[0:80])
                #return(True)
        else:
            #raise Exception("Fasta file can not be read.",self.filepath)
            print("Fasta file can not be read.",self.filepath)
    def read_sequences(self):
        "Generate rows for insertion in database"
        for seq_record in SeqIO.parse(self.filepath, "fasta"):
            parts = seq_record.id.split('|')
            if len(parts) < 2:
                raise Exception('Missing metadata',seq_record.id, self.filepath)
            #node_id = self.pldb.get_locus_uid(parts[0])
            reversed=re.match('^_R_',parts[0]) is not None
            req_id = dict(map(lambda x: tuple(x.split('_')),parts[2:]))
            yield self.build_row(req_id,reversed,seq_record.seq)
    def upload_sequences(self):
        #raise NotImplementedError('Sequence upload not yet implented')
        if self.file_uid is None:
            self.load_file_ids()
        if self.pldb.conn is None:
            self.pldb.open_connection()
            #self.load_file_ids()
        #print(self.insRow);
        #raise NotImplementedError('Sequence upload not yet implented')
        #raise NotImplementedError('Sequence scan not yet implented')
        if os.stat(self.filepath)[6] > 3:
            c=self.pldb.conn.cursor(prepared=True)
            with self.pldb.ww.lock:
                if not self.pldb.ww.ticking:
                    self.pldb.ww.resume_place()
            tick_cycle = itertools.repeat(True) # ( i % 10 == 0 for i in itertools.count(1) ) )
            for row_set in itertools.izip_longest(*(self.read_sequences(),)*10):
                if not self.pldb.dry_run:
                    c.executemany(self.insRow,list(itertools.takewhile(lambda x: x is not None, row_set)) )
                    if next( tick_cycle):
                        if c.rowcount >0:
                            print(10-c.rowcount,end='',file=self.pldb.ww)
                        else:
                            print('.',end='',file=self.pldb.ww)
                else:
                    if next( tick_cycle):
                        print('U',end='',file=self.pldb.ww)
                if self.aln_len is None and row_set[0] is not None:
                    self.aln_len = len(row_set[0][3])
            if self.aln_len is not None:
                if not self.pldb.dry_run:
                    c.execute(self.insLen,(self.file_uid,self.aln_len))
                    if c.rowcount >0:
                        print('#',end='',file=self.pldb.ww)
                    else:
                        print('%',end='',file=self.pldb.ww)
                else:
                    print('L',end='',file=self.pldb.ww)
            print('^',end='',file=sys.stderr)
            print('^',end='',file=self.pldb.ww)
            return([(self.filepath,1)])
        else:
            #raise Exception("Fasta file can not be read.",self.filepath)
            print("Fasta file can not be read.",self.filepath)
            return([])

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)

def wrap(bytes,n=70):
    g=grouper(bytes,70,'')
    return itertools.imap(lambda x: ''.join(x),g)


