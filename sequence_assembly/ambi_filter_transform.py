#!/usr/bin/env python2
##useful

from __future__ import print_function
from ambi_sites_filter import AbstractFilter
from tio import AlignmentTransform,TyranniSrc,get_wedge_args
from sites import SiteSeqCounter
import tyranniSQL

class RemoteFilter(AbstractFilter,AlignmentTransform):
    def __init__(self,verbose=False,pd_query=None,**kwargs):
        AbstractFilter.__init__(self,verbose=verbose)
        AlignmentTransform.__init__(self,**kwargs)
        if pd_query is not None:
            self.pd_q = pd_query
        else:
            self.pd_q = tyranniSQL.get_connection_from_args(tyranniSQL.DEFAULT)
    def query(self):
        return self.pd_q.named_cursor()
    scan_private_alleles="""SELECT * from private_allelles where aln_uid=%(aln_uid)s"""
    count_private_alleles="""SELECT count(*) n from private_allelles 
    where aln_uid=%(aln_uid)s"""
    def has_private_alleles(self,aln_uid,min_count=4):
        c=self.query()
        c.execute(self.count_private_alleles,{'aln_uid':aln_uid})
        res=c.fetchall()
        return res[0].n >= min_count
    def get_private_alleles(self,aln_uid):
        c=self.query()
        c.execute(self.scan_private_alleles,{'aln_uid':aln_uid})
        return sorted([r.site_num for r in iter(c)])


# def get_filter_args(parser):
#     parser.add_argument(
#             "-s","--slice",metavar="num_slices",default=8,
#             type=int,
#             help="""Process files from one of num_slices sets [%(default)s]""")
#     parser.add_argument(
#             "-w","--wedge",type=int, default=0,
#             help="""The file set from which to draw [%(default)s]""")
#     parser.add_argument(
#             "-f","--src-file",type=int,
#             help="""A single file to operate on""")
#     return parser
#     

def get_args():
    import argparse
    parser = argparse.ArgumentParser(
            description="""Filter alignments in batchgroup""")
    parser.add_argument(
            "-v","--verbose",default=False,
            action="store_true",
            help="""Output diagnostics""")
    w = parser.add_argument_group('Wedge Selection')
    get_wedge_args(w)
    tmp = parser.add_mutually_exclusive_group(required=True)
    
    tmp.add_argument(
            "--list-groups",default=False,
            action="store_true",
            help="""List known batch groups""")
    tmp.add_argument(
            "-i","--input-group",type=str,
            help="""The batch group from which to draw""")
    parser.add_argument(
            "-o","--output-group",default="NEW_GROUP",type=str,
            help="""The batch group of the result [%(default)s]""")
    parser.add_argument(
            "--cores",default=4,type=int,
            help="""The number of workers to use [%(default)s]""")
    return parser.parse_args()


def get_feeder_args(io_args):
    import argparse
    return(argparse.Namespace(src_col='file_uid',
    wedge_col='locus_uid',
    query="""select aln_file_uid file_uid, locus_uid
FROM stage5jobs join jobs using(job_uid) join alignment_length on aln_file_uid = file_uid
where batchgroup='{0.input_group}' and len > 0""".format(io_args)))

#join alignment_depth_cache on aln_file_uid=file_uid
# class BatchGroupSrc(TyranniSrc):
#     def __init__(self,**kwargs):
#         TyranniSrc.__init__(self,**kwargs)
# #     def get_cmd(self):
# #         if self.args is None:
# #             return(process_src_cmd(self.cmd))
# #         elif self.args.src_file is not None:
# #             wedge_cmd="""{1.src_col}={0.src_file} """.format(self.args,self.cmd)
# #             return(process_src_cmd(self.cmd,wedge_cmd))
# #         else:
# #             wedge_cmd="""mod({1.wedge_col},{0.slice})={0.wedge} 
# # and mod(floor({1.wedge_col}/{0.slice}),{0.count})=%s """.format(self.args,self.cmd)
# #             return(process_src_cmd(self.cmd,wedge_cmd))
#         

def worker(src, args, ticker=None):
    import sys
    filter_worker=RemoteFilter(group_args=args, ticker=ticker)
    stat_worker=SiteSeqCounter(pd_in=filter_worker.pd_out)
    #pd_q = tyranniSQL.get_connection_from_args(tyranniSQL.DEFAULT)
    #c=pd_q.named_cursor()
    for f in iter(src.get, 'STOP'):
        try:
            tmp=filter_worker.trim_file_no_stats(f)
            filter_worker.commit()
            stat_worker.scan_file(tmp)
            stat_worker.commit()
            print('S',end='',file=ticker)
        except Exception:
            print('!',end='',file=ticker)
            pass
        sys.stdout.flush()
    src.close()
    filter_worker.close()
    stat_worker.close()
    src.join_thread()

def list_groups():
    pd_q = tyranniSQL.get_connection_from_args(tyranniSQL.DEFAULT)
    c=pd_q.named_cursor()
    c.execute("""SELECT distinct batchgroup 
    from stage5jobs join jobs using(job_uid)""")
    for r in iter(c):
        print(r)
    return

def mp_main():
    args=get_args()
    if args.list_groups:
        list_groups()
        return
    if args.slice>1:
        src=TyranniSrc(get_feeder_args(args),wedge_args=args)
    else:
        src=TyranniSrc(get_feeder_args(args))
    # src.pd_in.ww.wrap=30
    ticker=src.pd_in.ww
    ticker.resume_place()
    job_exists="""SELECT aln_file_uid file_uid from  stage5jobs join jobs using(job_uid)
where batchgroup='{0.output_group}' and locus_uid = %(locus_uid)s""".format(args)
    files=src.get_rows()
    from multiprocessing import Process,Queue
    jobs=Queue()
    # w0=Process(target=worker,args=(jobs,args))
    # w1=Process(target=worker,args=(jobs,args))
    # w2=Process(target=worker,args=(jobs,args))
    workers = [ Process(target=worker,args=(jobs,args,ticker)) for _ in range(args.cores)]
    # w0.start(); w1.start(); w2.start()
    for w in workers:
        w.start()
    pd_q = tyranniSQL.get_connection_from_args(tyranniSQL.DEFAULT,ticker=ticker)
    c=pd_q.named_cursor()
    for f in files:
        c.execute(job_exists,f._asdict())
        if len(c.fetchall()) >0:
            print("Found output for", f)
            continue
        jobs.put(f.file_uid)
    for _ in workers:
        jobs.put('STOP')
    # jobs.put('STOP')
    # jobs.put('STOP')
    jobs.close()
    jobs.join_thread()
    for w in workers:
        w.join()
    ticker.break_count()

def main():
    args=get_args()
    if args.list_groups:
        list_groups()
        return
    filter_worker=RemoteFilter(group_args=args)
    stat_worker=SiteSeqCounter(pd_in=filter_worker.pd_out)
    if args.slice>1:
        src=TyranniSrc(get_feeder_args(args),wedge_args=args)
    else:
        src=TyranniSrc(get_feeder_args(args))
    job_exists="""SELECT aln_file_uid file_uid from  stage5jobs join jobs using(job_uid)
where batchgroup='{0.output_group}' and locus_uid = %(locus_uid)s""".format(args)
    files=src.get_rows()
#    f=next(files)
#    print(f)
#    tmp=filter_worker.trim_file_no_stats(f)
#    stat_worker.scan_file(tmp)
    pd_q = tyranniSQL.get_connection_from_args(tyranniSQL.DEFAULT)
    c=pd_q.named_cursor()
    for f in files:
        #print(f)
        c.execute(job_exists,f._asdict())
        if len(c.fetchall()) >0:
            print("Found output for", f)
            continue
        tmp=filter_worker.trim_file_no_stats(f.file_uid)
        #print(tmp)
        #filter_worker.pd_in.ww.flush()
        filter_worker.commit()
        try:
            stat_worker.scan_file(tmp)
            stat_worker.commit()
        except Exception:
            print('!',end='',file=src.pd_in.ww)
            pass

if __name__ == '__main__':
    mp_main()
