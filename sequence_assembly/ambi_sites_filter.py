#!/usr/bin/env python2
##useful

from __future__ import print_function
import re
from Bio import SeqIO,Seq,SeqRecord
from collections import Counter
ambi=re.compile('[MRWSKY][^MRWSKY]*[MRWSKY]([^MRWSKY]*[MRWSKY])+')
nongap=re.compile('[^-?]')
acgt20=re.compile('([^ACGTMRWSKY]*[ACGT]){20,}')



def main():
    import argparse
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "-v","--verbose",default=False,
            action="store_true",
            help="""Output diagnostics"""
        )
    parser.add_argument(
            "infile",
            type=str,
            help="""The input file"""
        )
    parser.add_argument(
            "outfile",
            type=str,
            help="""The output file"""
        )
    args=parser.parse_args()
    worker=LocalFilter(args.verbose)
    write_to_file(worker.get_seqs(args.infile),args.outfile)


class AbstractFilter(object):
    def __init__(self,verbose=False):
        self._v=verbose
        self.delta=0
        self.changed=0
        self._name=None
    def has_private_alleles(self,aln_uid,min_count=4):
        raise NotImplementedError()
    def get_private_alleles(self,aln_uid):
        raise NotImplementedError()
    def reset_counts(self):
        self.delta=0
        self.changed=0
    def trim_seq(self,tmp_rec):
        src_seq=tmp_rec.aligned_seq_data
        foo=ambi.search(str(src_seq))
        if foo is not None:
#             print(tmp_rec.aln_uid, "Removing ambiguous")
            seq=sub_ambi(foo,edit_seq=src_seq)
        else:
            seq=src_seq
        if self.has_private_alleles(tmp_rec.aln_uid,4):
#             print(tmp_rec.aln_uid, "Scan Private")
            seq,drop=sub_priv(src_seq,
                                 self.get_private_alleles(tmp_rec.aln_uid),seq)#,
#                        tmp_rec.aligned_seq_data.alphabet)
            if len(drop) >0:
                if min(drop)>1 and get_acgt_count(src_seq[0:min(drop)-1])<10:
                    seq=subst_N_sub(seq,endpos=min(drop)-1,mutable=True,subst='-')
                if max(drop)<len(seq) and get_acgt_count(src_seq[max(drop):])<10:
                    seq=subst_N_sub(seq,pos=max(drop),mutable=True,subst='-')
        tmp_delta=0
        if seq is not src_seq:
            for a,b in zip(src_seq,seq):
                if a != b:
                    self.delta+=1
                    tmp_delta+=1
                    #print('N',end='')
        if tmp_delta > 0:
            if self._v: 
                print('N'*tmp_delta)
            self.changed+=1
            return tmp_rec._replace(aligned_seq_data=seq)
        return tmp_rec
    def trim_all(self, seqs):
        for s in seqs:
            yield self.trim_seq(s)
        if self._v:
            print(self._name,"delta",self.delta,"bp")
            print(self._name,"changed",self.changed,"sequences")


class LocalFilter(AbstractFilter):
    def __init__(self,verbose=False):
        AbstractFilter.__init__(self,verbose=verbose)
        from sites_local import SiteCounter
        self.proc=SiteCounter()
    def has_private_alleles(self,aln_uid,min_count=4):
        return self.proc.has_private_alleles(aln_uid,min_count=min_count)
    def get_private_alleles(self,aln_uid):
        return self.proc.get_private_alleles(aln_uid)
    def get_seqs(self, file_path,file_type="fasta"):
        from sites_local import read_sequences_alt
        self.proc.scan_file(file_path,file_type)
        if self._v: 
            print("Sites Scanned")
            import os.path
            self._name=os.path.basename(file_path)
        self.reset_counts()
        return self.trim_all(read_sequences_alt(file_path,file_type))

def write_to_file(recs, file_path,file_type="fasta",strip_aln_uid=True):
    from itertools import imap
    if strip_aln_uid:
        seqs=imap(lambda r: SeqRecord.SeqRecord(id='|'.join([r.aln_uid]+[ p for p in r.s_info if not p.startswith('ALNUID')]),
                                                description='',seq=r.aligned_seq_data),
                  recs)
    else:
        seqs=imap(lambda r: SeqRecord.SeqRecord(id='|'.join([r.aln_uid]+r.s_info),description='',
                                                seq=r.aligned_seq_data),
                  recs)
    if file_type != "fasta":
        SeqIO.write(seqs,file_path,file_type)
        return
    from itertools import izip_longest
    with open(file_path,'wb') as h:
        for s in seqs:
            print('>{}'.format(s.id),file=h)
            c = iter(s.seq)
            for l in izip_longest(*[c]*70,fillvalue=''):
                print(*l,sep='',file=h)


def test_main2():
    fList=["aln-20161112PASS.old/gblocks_uce-11.fasta","aln-20161112PASS.old/gblocks_uce-11.fasta_out2"]
    outList=["trim_gblocks_uce-11.fasta","trim_gblocks_uce-11.fasta_out2"]
    worker=LocalFilter(True)
    write_to_file(worker.get_seqs(fList[0]),outList[0])
    write_to_file(worker.get_seqs(fList[1]),outList[1])

def test_main():
    fList=["aln-20161112PASS.old/gblocks_uce-11.fasta","aln-20161112PASS.old/gblocks_uce-11.fasta_out2"]
    # for testing
    from sites_local import SiteCounter, read_sequences_alt
    proc=SiteCounter()
    proc.scan_file(fList[0],"fasta")
    print("Sites Scanned")
    import sys
    sys.stdout.flush()
    delta=0
    changed=0
    for tmp_rec in read_sequences_alt(fList[0],"fasta"):
        src_seq=tmp_rec.aligned_seq_data
        foo=ambi.search(str(src_seq))
        if foo is not None:
#             print(tmp_rec.aln_uid, "Removing ambiguous")
            seq=Seq.Seq(sub_ambi(foo),src_seq.alphabet)
        else:
            seq=src_seq
        if proc.has_private_alleles(tmp_rec.aln_uid,4):
#             print(tmp_rec.aln_uid, "Scan Private")
            seq,drop=sub_priv(src_seq,
                                 proc.get_private_alleles(tmp_rec.aln_uid),seq)#,
#                        tmp_rec.aligned_seq_data.alphabet)
            if len(drop) >0:
                if min(drop)>1 and get_acgt_count(src_seq[0:min(drop)-1])<10:
                    seq=subst_N_sub(seq,endpos=min(drop)-1,mutable=True)
                if max(drop)<len(seq) and get_acgt_count(src_seq[max(drop):])<10:
                    seq=subst_N_sub(seq,pos=max(drop),mutable=True)
        tmp_delta=0
        if seq is not src_seq:
            for a,b in zip(src_seq,seq):
                if a != b:
                    delta+=1
                    tmp_delta+=1
                    #print('N',end='')
        if tmp_delta > 0:
            print('N'*tmp_delta)
            sys.stdout.flush()
            changed+=1
#         else:
#             print(tmp_rec.aln_uid, "Unchanged")
    print(fList[0],"delta",delta,"bp")
    print(fList[0],"changed",changed,"sequences")
    proc.scan_file(fList[1],"fasta")
    print("Sites Scanned")
    sys.stdout.flush()
    delta=0
    changed=0
    for tmp_rec in read_sequences_alt(fList[1],"fasta"):
        foo=ambi.search(str(tmp_rec.aligned_seq_data))
        if foo is not None:
#             print(tmp_rec.aln_uid, "Removing ambiguous")
            seq=Seq.Seq(sub_ambi(foo),tmp_rec.aligned_seq_data.alphabet)
        else:
            seq=tmp_rec.aligned_seq_data
        if proc.has_private_alleles(tmp_rec.aln_uid,4):
#             print(tmp_rec.aln_uid, "Scan Private")
            seq,drop=sub_priv(tmp_rec.aligned_seq_data,
                                 proc.get_private_alleles(tmp_rec.aln_uid),seq)#,
#                        tmp_rec.aligned_seq_data.alphabet)
        tmp_delta=0
        if seq is not tmp_rec.aligned_seq_data:
            for a,b in zip(tmp_rec.aligned_seq_data,seq):
                if a != b:
                    delta+=1
                    tmp_delta+=1
                    #print('N',end='')
        if tmp_delta > 0:
            print('N'*tmp_delta)
            sys.stdout.flush()
            changed+=1
#         else:
#             print(tmp_rec.aln_uid, "Unchanged")
    print(fList[1],"delta",delta,"bp")
    print(fList[1],"changed",changed,"sequences")
#     for tmp_rec in SeqIO.parse(fList[1],"fasta"):
#         foo=ambi.search(str(tmp_rec.seq))
#         if foo is not None:
#             print(tmp_rec.id)
#             print(sub_ambi(foo))


def get_counts(rec_seq):
    bar=dict(Counter(str(rec_seq)))
    return((sum((v for k,v in bar.iteritems() if k in 'ACGT')),
            sum((v for k,v in bar.iteritems() if k in 'MRWSKY'))))

def get_acgt_count(rec_seq):
    return(sum((v for k,v in Counter(str(rec_seq)).iteritems() if k in 'ACGT')))

def subst_N_all(seq, subst='N'):
    if not isinstance(seq,str):
        seq=str(seq)
    return nongap.sub(subst, seq)

def subst_N_sub(seq, pos=None, endpos=None, mutable=False, subst='N'):
    if mutable:
        if seq.__class__=="Bio.Seq.MutableSeq":
            pass
        elif seq.__class__=="Bio.Seq.Seq":
            import Bio.Seq.MutableSeq
            seq=Bio.Seq.MutableSeq(seq,seq.alphabet)
        elif isinstance(seq,bytearray):
            pass
        else:
            mutable=False
    if pos is not None:
        if endpos is not None:
            subseq=seq[pos:endpos]
            if mutable:
                seq[pos:endpos]=nongap.sub(subst,str(subseq))
                return(seq)
        else:
            subseq=seq[pos:]
            if mutable:
                seq[pos:]=nongap.sub(subst,str(subseq))
                return(seq)
    else:
        if endpos is None:
            if mutable:
                seq[:]=nongap.sub(subst,str(subseq))
                return(seq)
            return subst_N_all(seq, subst=subst)
        subseq=seq[:endpos]
    if not isinstance(seq,str):
        subseq=str(subseq)
    res=nongap.sub(subst, subseq)
    if pos is not None:
        res = seq[:pos]+res
    if endpos is not None:
        res = res + seq[endpos:]
    return res

def sub_ambi(foo,l_sect='',r_sect='',subst='N',edit_seq=None):
    if edit_seq is None:
        edit_seq=foo.string
    good,poor=get_counts(foo.group(0))
    if good <= 19* poor:
        #nongap.sub('N',foo.group(0))
        if l_sect is not None:
            l_sect+=edit_seq[0:foo.start()]
            if get_acgt_count(l_sect)<10:
                #l_sect=nongap.sub(subst,str(l_sect))
                l_sect=subst_N_all(l_sect,subst)
        else:
            l_sect=edit_seq[0:foo.start()]
        if r_sect is not None:
            r_sect=edit_seq[foo.end():]+r_sect
            if get_acgt_count(r_sect)<10:
                r_sect=subst_N_all(r_sect,subst)
                #r_sect=nongap.sub(subst,str(r_sect))
        else:
            r_sect=edit_seq[foo.end():]
        return(l_sect+nongap.sub(subst,foo.group(0))+r_sect)
    elif poor<4:
        return(edit_seq)
    parM=foo.group(0)
    middle=acgt20.search(parM)
    if middle is None:
        return(edit_seq)
    res=list()
    parL=edit_seq[0:foo.start()]
    seqL=parM[0:middle.start()]
    l_srch=ambi.search(str(seqL))
    if l_srch is None:
        if l_sect is not None:
            res.append(l_sect)
        res.append(parL)
        res.append(seqL)
    elif l_sect is None:
        res.append(parL)
        res.append(sub_ambi(l_srch,None,None,subst=subst,edit_seq=seqL))
    else:
        l_sect+=parL
        if get_acgt_count(l_sect)<10:
            res.append(sub_ambi(l_srch,l_sect,None,subst=subst,edit_seq=seqL))
        else:
            res.append(l_sect)
            res.append(sub_ambi(l_srch,None,None,subst=subst,edit_seq=seqL))
    res.append(middle.group(0))
    parR=edit_seq[foo.end():]
    seqR=parM[middle.end():]
    r_srch=ambi.search(str(seqR))
    if r_srch is None:
        res.append(seqR)
        res.append(parR)
        if r_sect is not None:
            res.append(r_sect)
    elif r_sect is None:
        res.append(sub_ambi(r_srch,None,None,subst=subst,edit_seq=seqR))
        res.append(parR)
    else:
        r_sect=parR+r_sect
        if get_acgt_count(r_sect)<10:
            res.append(sub_ambi(r_srch,None,r_sect,subst=subst,edit_seq=seqR))
        else:
            res.append(sub_ambi(r_srch,None,None,subst=subst,edit_seq=seqR))
            res.append(r_sect)
    return(concat(res))

def concat(objs):
    res=None
    for i in objs:
        if res is None:
            res=i
        else:
            res=res+i
    return res

def sub_priv(seq,priv_sites,edit_seq=None, subst='-'):
    if edit_seq is None:
        edit_seq=seq
    if len(priv_sites)<4:
        return(edit_seq,list())
    priv_max=max(priv_sites)
    priv_min=min(priv_sites)-1
    if 4*len(priv_sites)>=get_acgt_count(seq[priv_min:priv_max]):
        return((edit_seq[0:priv_min]+
                        nongap.sub(subst, str(edit_seq[priv_min:priv_max]))+
                        edit_seq[priv_max:],priv_sites))
    #keep_sites=list()
    cur_sites=list()
    segs=list()
    cur_l =cur_r =0;
    while len(priv_sites)>3:
        while get_acgt_count(seq[(priv_sites[0]-1):priv_sites[3]])>16:
            priv_sites.pop(0)
            if len(priv_sites)<4:
                if cur_r == 0:
                    return((edit_seq,list()))
                seg_iter=iter(segs)
                seg=next(seg_iter)
                for i in seg_iter:
                    seg+=i
                return((seg+edit_seq[cur_r:],cur_sites))
        cur_l=priv_sites[0]-1
        segs.append(seq[cur_r:cur_l])
        for _ in range(4):
            cur_sites.append(priv_sites.pop(0))
        if len(priv_sites)>0:
            while get_acgt_count(seq[cur_l:priv_sites[0]])<=4+4*len(cur_sites):
                cur_sites.append(priv_sites.pop(0))
                if len(priv_sites)<1:
                    break
        cur_r=cur_sites[-1]
        segs.append(nongap.sub(subst, str(edit_seq[cur_l:cur_r])))
    seg_iter=iter(segs)
    seg=next(seg_iter)
    for i in seg_iter:
        seg+=i
    return(seg+edit_seq[cur_r:],cur_sites)


if __name__ == '__main__':
    main()
