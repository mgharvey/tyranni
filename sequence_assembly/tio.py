#!/usr/bin/env python
from __future__ import print_function

import tyranniSQL
import itertools,re
from collections import namedtuple
import multiprocessing
import tyranniSQL.helpers
# from argparse import Namespace

class TyranniIO(object):
    def __init__(self, pd_in=None, pd_out=None, ticker=None):
        if pd_in is not None:
            self.pd_in = pd_in
        else:
            self.pd_in = tyranniSQL.get_connection_from_args(tyranniSQL.DEFAULT,ticker=ticker)
        if pd_out is not None:
            self.pd_out = pd_out
        else:
            self.pd_out = tyranniSQL.get_connection_from_args(tyranniSQL.COMMIT,ticker=self.pd_in.ww)
    def reader(self):
        return self.pd_in.named_cursor()
    def writer(self):
        return self.pd_out.named_cursor()
    def writer_prepared(self):
        return self.pd_out.conn.cursor(prepared=True)
    def commit(self):
        self.pd_in.conn.commit()
        return self.pd_out.conn.commit()
    def rollback(self):
        self.pd_in.conn.rollback()
        return self.pd_out.conn.rollback()
    def close(self):
        self.pd_in.conn.close()
        return self.pd_out.conn.close()
    def rolling_command(self, cmd, args,
                        cur=None, named=True,
                        line_commit=False, term_commit=False):
        if cur is None:
            cur = self.writer()
        if named:
            if line_commit:
                for arg in args:
                    cur.execute(cmd,arg._asdict())
                    yield arg
                    self.commit()
            else:
                for arg in args:
                    cur.execute(cmd,arg._asdict())
                    yield arg
        else:
            if line_commit:
                for arg in args:
                    cur.execute(cmd,arg)
                    yield arg
                    self.commit()
            else:
                for arg in args:
                    cur.execute(cmd,arg)
                    yield arg
        if term_commit:
            self.commit()

def get_wedge_args(parser):
#     parser = argparse.ArgumentParser(
#             description="""Trim alignments for a set of files""")
    parser.add_argument(
            "-c","--count",metavar="sub_slices",default=8,
            type=int,
            help="""The number of sub slices [%(default)s]""")
    parser.add_argument(
            "-s","--slice",metavar="num_slices",default=8,
            type=int,
            help="""Process files from one of num_slices sets [%(default)s]""")
    parser.add_argument(
            "-w","--wedge",type=int, default=0,
            help="""The file set from which to draw [%(default)s]""")
    parser.add_argument(
            "-f","--src-file",type=int,
            help="""A single file to operate on""")
    return parser#.parse_args()




def process_src_cmd(cmd_args,wedge_cmd=None):
    j_cmd="""SELECT {0.src_col} from {0.src_table} left join 
{0.sink_table} a using ({0.joint_col}) 
where a.{0.joint_col} is NULL"""
    s_cmd="""SELECT {0.src_col} from {0.src_table} left join 
{0.sink_table} a using ({0.src_col}) 
where a.{0.joint_col} is NULL"""
    l_cmd="limit {0.limit_rows}"
    if 'query' in dir(cmd_args):
        cmd=cmd_args.query
    elif 'joint_col' in dir(cmd_args):
        cmd=j_cmd
    else:
        cmd=s_cmd
    if wedge_cmd is not None:
        cmd=' and '.join((cmd,wedge_cmd))
    if 'limit_rows' in dir(cmd_args):
        cmd=' '.join((cmd,l_cmd))
    return(cmd.format(cmd_args))


class TyranniSrc(object):
    def __init__(self, base_cmd, auto_tick='+', wedge_args=None, pd_in=None, ticker=None):
        if pd_in is not None:
            self.pd_in = pd_in
        else:
            self.pd_in = tyranniSQL.get_connection_from_args(tyranniSQL.DEFAULT,ticker=ticker)
        self.auto_tick=auto_tick
        self.cmd=base_cmd
        self.args=wedge_args
    def reader(self):
        return self.pd_in.named_cursor()
    def close(self):
        self.pd_in.ww.break_count()
        self.pd_in.conn.close()
    def tick(self):
        if self.auto_tick is not None:
            print(self.auto_tick,end='',file=self.pd_in.ww)
    def get_cmd(self):
        if self.args is None:
            return(process_src_cmd(self.cmd))
        elif self.args.src_file is not None:
            wedge_cmd="""{1.src_col}={0.src_file} """.format(self.args,self.cmd)
            return(process_src_cmd(self.cmd,wedge_cmd))
        else:
            wedge_cmd="""mod({1.wedge_col},{0.slice})={0.wedge} 
and mod(floor({1.wedge_col}/{0.slice}),{0.count})=%s """.format(self.args,self.cmd)
            return(process_src_cmd(self.cmd,wedge_cmd))
    def get_data(self):
        sql_cmd=self.get_cmd()
        if sql_cmd is None:
            yield self.args.src_file
            self.tick()
            return
        if self.args is None:
            c=self.reader()
            c.execute(sql_cmd)
            for r in iter(c):
                yield r[0]
                self.tick()
            return
        c=self.reader()
        for sub_slice in range(self.args.count):
            c.execute(sql_cmd,(sub_slice,))
            for r in iter(c):
                yield r[0]
                self.tick()
        return
    def get_rows(self):
        sql_cmd=self.get_cmd()
        if sql_cmd is None:
            yield (self.args.src_file,)
            self.tick()
            return
        if self.args is None:
            c=self.reader()
            c.execute(sql_cmd)
            for r in iter(c):
                yield r
                self.tick()
            return
        c=self.reader()
        for sub_slice in range(self.args.count):
            c.execute(sql_cmd,(sub_slice,))
            for r in iter(c):
                yield r
                self.tick()
        return

class AlignmentTransform(TyranniIO):
    def __init__(self,group_args,**kwargs):
        TyranniIO.__init__(self,**kwargs)
        self.io_args=group_args
    site_scanner_lim="""SELECT site_num,snp_count,ACGT_total FROM alignment_depth_cache 
where file_uid=%s and ACGT_total>=%s limit %s"""
    site_scanner="""SELECT site_num,snp_count,ACGT_total FROM alignment_depth_cache 
where file_uid=%s and ACGT_total>=%s"""
    site_scanner_full_lim="""SELECT * FROM alignment_depth_cache 
where file_uid=%s and ACGT_total>=%s limit %s"""
    site_scanner_full="""SELECT * FROM alignment_depth_cache 
where file_uid=%s and ACGT_total>=%s"""
    def scan_site_counts(self,in_file,lim=20,cur=None):
        if cur is None:
            c=self.reader()
        else:
            c=cur
        if lim is not None:
            c.execute(self.site_scanner_lim,(in_file,self.io_args.min_depth,lim))
        else:
            c.execute(self.site_scanner,(in_file,self.io_args.min_depth))
        return(iter(c))
    def scan_sequences(self,in_file,lim=100,cur=None):
        if cur is None:
            c=self.reader()
        else:
            c=cur
        if lim is None:
            c.execute(self.seq_scanner,(in_file,))
        else:
            c.execute(self.seq_scanner_lim,(in_file,lim))
        return(iter(c))
    seq_scanner_lim="""SELECT * FROM aligned_sequences where file_uid=%s limit %s"""
    seq_scanner="""SELECT * FROM aligned_sequences where file_uid=%s """
    job_entry="""INSERT INTO jobs (ghost_ref, job_name, stage, Branch, Batch, BatchGroup, job_status) VALUES(%(ghost_ref)s, %(job_name)s, %(stage)s, %(Branch)s, %(Batch)s, %(BatchGroup)s, %(job_status)s)"""
    file_entry="""INSERT INTO files (file_category, file_hash, file_status, file_size, job_uid_temp, file_name_temp) VALUES(%(file_category)s, %(file_hash)s, %(file_status)s, %(file_size)s, %(job_uid)s, %(file_name)s)"""
    path_entry="""INSERT INTO job_files (file_uid, job_uid, file_name, purpose) VALUES(%(file_uid)s, %(job_uid)s, %(file_name)s, %(purpose)s)"""
    s5job_entry="""INSERT INTO stage5jobs (job_uid, eff_job_id, locus_uid, seq_file_uid, aln_file_uid, old_aln_file_uid) VALUES(%(job_uid)s, %(job_uid)s,%(locus_uid)s, %(seq_file_uid)s, %(file_uid)s, %(old_aln_file_uid)s)"""
    full_scanner="""SELECT s.locus_uid, j.job_name,j.Batch,p.file_name
FROM stage5jobs s join jobs j using(job_uid) 
join job_files p on p.job_uid=j.job_uid and  p.file_uid=s.aln_file_uid 
where aln_file_uid=%s and BatchGroup=%s limit 2"""
    aln_entry="""INSERT INTO aligned_sequences (seq_uid, sample_uid, locus_uid, file_uid, reversed)
SELECT seq_uid, sample_uid, locus_uid, %s, reversed from aligned_sequences where file_uid=%s"""
    def insert_file_entry(self,in_file_uid):
        out_data={'ghost_ref':None,'seq_file_uid':None,'stage':u'stage5',
                  'Branch':u'normal','job_status':0,'purpose':1,
                  'file_category':u'aligned_sequences', 'file_hash':None,
                  'file_status':0, 'file_size':None}
        cR=self.reader()
        cR.execute(self.full_scanner,(in_file_uid,self.io_args.input_group))
        rows=cR.fetchall()
        if cR.rowcount >1:
            print("Warning, multiple entries detected:",rows[0],rows[1])
        if cR.rowcount ==0:
            raise IndexError("Warning, no entries detected:",self.io_args,in_file_uid)
        in_job_full=rows[0]
        out_data['BatchGroup']=self.io_args.output_group
        out_data['locus_uid']=in_job_full.locus_uid
        out_data['old_aln_file_uid']=in_file_uid
        out_data['job_name']=in_job_full.job_name
        out_data['Batch']=re.sub(self.io_args.input_group,self.io_args.output_group,in_job_full.Batch)
        out_data['file_name']=in_job_full.file_name
        out_data['job_uid']=None
        out_data['file_uid']=None
        cW=self.writer()
        cW.execute(self.job_entry,out_data)
        out_data['job_uid']=cW.lastrowid
        cW.execute(self.file_entry,out_data)
        out_data['file_uid']=cW.lastrowid
        cW.execute(self.path_entry,out_data)
        cW.execute(self.s5job_entry,out_data)
        cW.execute(self.aln_entry,(out_data['file_uid'],in_file_uid))
        
        return( (out_data['file_uid'],out_data['job_uid']) )
    
    sql_norm="""UPDATE aligned_sequences SET aligned_seq_data=%(aligned_seq_data)s 
    WHERE seq_uid=%(seq_uid)s and file_uid=%(file_uid)s"""
    def upload_aligned_seq_data(self, rows, new_file_uid=None,
                                line_commit=False, term_commit=True):
        rows=iter(rows)
        if new_file_uid is not None:
            rows=itertools.imap(lambda r: r._replace(file_uid=new_file_uid), rows)
        res=self.rolling_command(self.sql_norm,rows,
                             line_commit=line_commit,
                             term_commit=term_commit)
        import mysql.connector.errors
        try:
            r=next(res)
            out_file=r.file_uid
            new_length=len(r.aligned_seq_data)
            count=1
            #self.pd_in.ww.write('U')
            for _ in res:
                #self.pd_in.ww.write('u')
                count+=1
            cOut=self.writer()
            cOut.execute("""INSERT INTO alignment_length (file_uid,len) 
                         VALUES (%s,%s)""",(out_file,new_length))
            #self.pd_in.ww.write('\n')
            return(count)
        except mysql.connector.errors.Error as e:
            print(e)
            print(e.stracktrace)
            return(0)
    def trim_seq(self,tmp_rec):
        raise NotImplementedError("Abstraction Trimmer!")
    def trim_all(self, seqs):
        for s in seqs:
            yield self.trim_seq(s)
    def trim_file_no_stats(self, src_file_uid):
        new_file_uid, new_job_uid=self.insert_file_entry(src_file_uid)
        src_rows=self.scan_sequences(src_file_uid,lim=None)
        row_count=self.upload_aligned_seq_data(self.trim_all(src_rows),
                                               new_file_uid=new_file_uid)
        #print(row_count)
        print('U',end='',file=self.pd_in.ww)
        return(new_file_uid)


    
