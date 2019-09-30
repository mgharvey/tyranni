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

Name: stage4.py

Author: Graham E Derryberry
Date: 11 September 2016

Description:
Compile stage 4 fasta files from the mysql db.


"""

class Stage4SamplesFile:
    """Represents a stage 4 fasta file to be downloaded from the database"""
    def __init__(self, pldb, locus_uid=None,locus_name=None,delay=False):
        if delay:
            raise NotImplementedException("Delay connection code not yet implemented")
        else:
            self.pldb=pldb
        self.pldb.ww.tick=50
        self.locus_name =locus_name
        self.locus_uid=locus_uid
        if self.locus_uid is not None:
            self.fetch_rows=''' select  sample_uid,locus_uid,seq_uid,branch,seq_data,depth_quality,sub_start,sub_end
            from best_unaligned_sequences join unaligned_sequences using (seq_uid,sample_uid,locus_uid)
            where locus_uid={}'''.format(self.locus_uid)
        else:
            self.fetch_rows=None
        if pldb.conn is not None:
            self.load_file_ids()
        else:
            self.file_ids=None
        # self.fetch_seq = ''' '''
    select_src_id='''select 
        file_uid, Order_no as sample_uid, sample_id as alias_uid, branch 
    from (select file_uid, job_uid, eff_job_id from known_output_files 
        where file_category='loci_for_sample' and file_path=%s ) as src
    join stage3stats as j on j.job_uid=src.job_uid or j.job_uid=src.eff_job_id'''
    def load_file_ids(self):
        "Get sample uid to sample name for stage 4 alignment dictionary"
        if self.pldb.conn is None:
            self.pldb.open_connection()
        self.file_ids=self.pldb.get_sample_dict(stage='stage4')
        if self.fetch_rows is None:
            self.fetch_rows=''' select  sample_uid,branch,locus_uid,seq_uid,seq_data,depth_quality,sub_start,sub_end
            from best_unaligned_sequences join unaligned_sequences using (seq_uid,sample_uid,locus_uid)
            where locus_uid={}'''.format(self.pldb.get_locus_uid(self.locus_name))
            self.locus_uid=self.pldb.get_locus_uid(self.locus_name)
    def download_sequences(self,out=sys.stdout,filter=None):
        if self.file_ids is None:
            if self.pldb.conn is None:
                self.pldb.open_connection()
            self.load_file_ids()
        
        c=self.pldb.named_cursor()
        #print(self.insRow);
        #raise NotImplementedError('Sequence upload not yet implented')
        #raise NotImplementedError('Sequence scan not yet implented')
        c.execute(self.fetch_rows)
        self.pldb.ww.resume_place()
        if filter is None:
            for row in c:
                print('>{}|{}|SUID_{}|LUID_{}|SEQUID_{}'.format(self.file_ids[row.sample_uid],row.branch,row.sample_uid,row.locus_uid,row.seq_uid),file=out)
                print(row.seq_data[row.sub_start:row.sub_end],file=out)
                print('+',end='',file=self.pldb.ww)
        else:
            for row in c:
                if not filter(row):
                    print('S',end='',file=self.pldb.ww)
                    continue
                print('>{}|{}|SUID_{}|LUID_{}|SEQUID_{}'.format(self.file_ids[row.sample_uid],row.branch,row.sample_uid,row.locus_uid,row.seq_uid),file=out)
                print(row.seq_data[row.sub_start:row.sub_end],file=out)
                print('+',end='',file=self.pldb.ww)



