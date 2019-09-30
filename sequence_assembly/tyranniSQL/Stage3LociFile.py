#!/usr/bin/env python2
##useful

from __future__ import print_function
import itertools, os
import datetime
from Bio import SeqIO
import mysql.connector
import re,os.path,sys,string
#import openpyxl,lxml
#from openpyxl import load_workbook
from mysql.connector.errors import IntegrityError
"""

Name: Stage3LociFile.py

Author: Graham E Derryberry
Date: 8 September 2016

Description:
Parse stage 3 fasta files for posting to mysql db.


"""

class Stage3LociFile:
    """Represents a stage 3 fasta file to be uploaded to the database"""
    def __init__(self, filepath,pldb):
        self.filepath=filepath
        self.pldb=pldb
        pldb.ww.tick=50
        if pldb.conn is not None:
            self.load_file_ids()
        else:
            self.file_ids=None
    select_src_id='''select 
        file_uid, Order_no as sample_uid, sample_id as alias_uid, branch 
    from (select file_uid, job_uid, eff_job_id from known_output_files 
        where file_category='loci_for_sample' and file_path=%s ) as src
    join stage3stats as j on j.job_uid=src.job_uid or j.job_uid=src.eff_job_id'''
    def load_file_ids(self):
        c=self.pldb.conn.cursor(named_tuple=True)
        c.execute(Stage3LociFile.select_src_id,(self.filepath,))
        res=c.fetchall()
        if c.rowcount >1:
            raise Exception('File identified multiple times',res, self.filepath)
        if c.rowcount <1:
            raise Exception('File not registered', self.filepath)
        self.file_ids=res[0]
        bonus=0
        if re.match('^normal$',self.file_ids.branch,flags=re.I) is not None:
            bonus=100
        elif re.match('^consensus$',self.file_ids.branch,flags=re.I) is not None:
            bonus=200
        self.insRow='''insert into unaligned_sequences 
        (sample_uid,alias_uid,locus_uid,file_uid,seq_data,branch,node,length,coverage,quality) values 
        ({},{},?,{},?,'{}',?,?,?,?+{:d})'''.format(
            self.file_ids.sample_uid,
            self.file_ids.alias_uid,
            self.file_ids.file_uid,
            self.file_ids.branch,
            bonus)
    reQ=re.compile('[^ACGTU]+')
    def assess_sequence(self,seq):
        return(2*len(Stage3LociFile.reQ.sub('',str(seq))))
    def build_row(self,locus_uid,seq,stats):
        return((locus_uid,
            format(seq),
            int(stats['NODE']),int(stats['length']),float(stats['cov']),
            2*len(Stage3LociFile.reQ.sub('',str(seq)))))
    def scan_sequences(self):
        if self.file_ids is None:
            if self.pldb.conn is None:
                self.pldb.open_connection()
            self.load_file_ids()
        print(self.insRow);
        #raise NotImplementedError('Sequence scan not yet implented')
        if os.stat(self.filepath)[6] > 3:
            for i,seq_record in itertools.izip(range(3),SeqIO.parse(self.filepath, "fasta")):
                parts = seq_record.id.split('|')
                if len(parts) < 2:
                    raise Exception('Missing metadata',seq_record.id, self.filepath)
                node_id = self.pldb.get_locus_uid(parts[0])
                stats = parts[1].split('_')
                s=dict(zip(stats[0::2],stats[1::2]))
                print(self.build_row(node_id,seq_record.seq,s))
                # print(dict(zip(stats[0::2],stats[1::2])))
                # print(seq_record.seq[0:80])
                #return(True)
        else:
            #raise Exception("Fasta file can not be read.",self.filepath)
            print("Fasta file can not be read.",self.filepath)
    def upload_sequences(self):
        if self.file_ids is None:
            if self.pldb.conn is None:
                self.pldb.open_connection()
            self.load_file_ids()
        print(self.insRow);
        #raise NotImplementedError('Sequence upload not yet implented')
        #raise NotImplementedError('Sequence scan not yet implented')
        if os.stat(self.filepath)[6] > 3:
            c=self.pldb.conn.cursor(prepared=True)
            c.execute('''select file_uid,count(*) from unaligned_sequences where file_uid=? group by file_uid''',
                (self.file_ids.file_uid,))
            tmp=c.fetchall()
            if c.rowcount >0:
                print(*tmp,sep='\n')
                return(tmp)
                # raise Exception("File already present",self.filepath) 
            #for i,seq_record in itertools.izip(range(3),SeqIO.parse(self.filepath, "fasta")):
            self.pldb.ww.resume_place()
            for seq_record in SeqIO.parse(self.filepath, "fasta"):
                parts = seq_record.id.split('|')
                if len(parts) < 2:
                    raise Exception('Missing metadata',seq_record.id, self.filepath)
                node_id = self.pldb.get_locus_uid(parts[0])
                stats = parts[1].split('_')
                s=dict(zip(stats[0::2],stats[1::2]))
                if not self.pldb.dry_run:
                    c.execute(self.insRow,self.build_row(node_id,seq_record.seq,s))
                # print(dict(zip(stats[0::2],stats[1::2])))
                # print(seq_record.seq[0:80])
                #return(True)
                print('+',end='',file=self.pldb.ww)
            c.execute('''select file_uid,count(*) from unaligned_sequences where file_uid=? group by file_uid''',
                (self.file_ids.file_uid,))
            res=c.fetchall()
            print(*res)
            return(res)
        else:
            #raise Exception("Fasta file can not be read.",self.filepath)
            print("Fasta file can not be read.",self.filepath)
            return([])


