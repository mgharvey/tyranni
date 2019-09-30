#!/usr/bin/env python2
##useful

from __future__ import print_function
import itertools, os
import datetime,random
import re,os.path,sys,string
from collections import namedtuple
from multiprocessing import RLock, Value
"""

Name: StdOutWordWrap.py

Author: Graham E Derryberry
Date: 11 September 2016

Description:
Basic class to visually handle single character output.

"""

PostBlock=namedtuple("PostBlock",["is_tick","lines"])
def empty_msg(is_tick=False):
    return PostBlock(is_tick=is_tick,lines=['']*is_tick)

class StdoutWordWrap(object):
    "Intended to mark progress with a wrapping column of characters"
    def __init__(self,wrap=100,tick=10):
        self.wrap = wrap
        self.tick = tick
        self.reNL=re.compile('\n')
        self.out=sys.stdout
        self.lock=RLock()
        self.ticking=Value('b',False,lock=self.lock)
        self.count = Value('i',0,lock=self.lock)
        self.l_cnt = Value('i',0,lock=self.lock)
    def item_count(self):
        "Get the number of progress ticks"
        with self.lock:
            return(self.count.value+self.l_cnt.value*self.wrap)
    def break_count(self):
        "Print total progress and reset counters"
        with self.lock:
            s='\nII: {:8d} items\n'.format(self.count.value+self.l_cnt.value*self.wrap)
            self.out.write(s)
            self.count.value= 0
            self.l_cnt.value = 0
            self.out.flush()
            self.ticking.value=False
            return(len(s))
    def resume_place(self):
        "Attempt to visually return to the correct column on screen"
        with self.lock:
            self.out.write('\n{:6d} '.format(self.l_cnt.value))
            self.out.write(' '*self.count.value)
            self.out.flush()
            self.ticking.value=True
            return(7+self.count.value)
    def write_msg(self,msg):
        with self.lock:
            if msg.is_tick:
                if not self.ticking:
                    self.resume_place()
                return(sum((self.write(l) for l in msg.lines )))
            if self.ticking.value:
                self.flush()
                self.out.write('\n')
                self.ticking.value=False
            iL=iter(msg.lines)
            self.out.write(next(iL))
            for l in iL:
                self.out.write('\n')
                self.out.write(l)
            self.flush()
            return(0)
    def write(self, stuff):
        """Act as a suitable target for print
        
        Each character in stuff increments the counter. 
        Each newline will reset the counter."""
        with self.lock:
            if self.reNL.match(stuff) is not None:
                lines = self.reNL.split(stuff,maxsplit=1)
                res=self.write(lines[0])+1
                self.out.write('\n')
                self.out.flush()
                self.count.value=0
                self.l_cnt.value=0
                if len(lines[1]) > 0:
                    res+=self.write(lines[1])
                else:
                    self.out.write('     0 ')
                    res+=7
                return(res)
            if len(stuff) + self.count.value > self.wrap:
                self.out.write(stuff[0:(self.wrap-self.count.value)])
                self.out.write('\n')
                self.l_cnt.value+=1
                self.out.write('{:6d} '.format(self.l_cnt.value))
                self.out.flush()
                res=self.wrap-self.count.value
                self.count.value=0
                res+=self.write(stuff[res:])+8
                return(res)
            self.out.write(stuff)
            if len(stuff) + self.count.value is self.wrap:
                self.out.write('\n')
                self.l_cnt.value+=1
                self.out.write('{:6d} '.format(self.l_cnt.value))
                self.out.flush()
                self.count.value=0
                return(len(stuff)+8)
            aTick = self.count.value/self.tick
            self.count.value+=len(stuff)
            if self.count.value/self.tick>aTick:
                self.out.flush()
            return(len(stuff))
    def flush(self):
        self.out.flush()

class ProgressTickerBuffer:
    "Buffer to reduce amount of locking done on StdoutWordWrap"
    def __init__(self,ww):
        self.ww=ww
        self.tick=random.choice([11,13,14,15,16,17,19])
        self.buf=''
    def write(self,stuff):
        "Write to buffer, post when sufficient data is loaded"
        self.buf+=stuff
        if len(self.buf) >=self.tick:
            self.ww.write(self.buf)
            self.buf=''
    def flush(self):
        self.ww.write(self.buf)
        self.buf=''
        self.ww.flush()


class MessageTickerBuffer:
    "Buffer to reduce amount of locking done on StdoutWordWrap"
    def __init__(self,ww):
        self.ww=ww
        self.first=empty_msg()
        self.last=empty_msg()
        self.msg_buf=list()
        self.count=0
        self.cur=None
        self.buf=None
    def write_msg(self,msg):
        if self.cur is not None:
            if self.cur.is_tick:
                self.count+=sum(map(len,self.buf))
            self.cur= None
            self.buf=None
        self.msg_buf.append(msg)
        if msg.is_tick:
            tmp=sum(map(len,msg.lines))
            self.count+=tmp
            return(tmp)
        return(0)
    def write(self,stuff):
        if self.cur is None:
            self.cur=empty_msg()
            self.buf=self.cur.lines
            self.msg_buf.append(self.cur)
        elif self.cur.is_tick:
            self.count+=sum(map(len,self.buf))
            self.cur=empty_msg()
            self.buf=self.cur.lines
            self.msg_buf.append(self.cur)
        self.buf.append(stuff)
        return(0)
    def tick(self,stuff):
        if self.cur is None:
            self.cur=empty_msg(True)
            self.buf=self.cur.lines
            self.msg_buf.append(self.cur)
        elif not self.cur.is_tick:
            self.cur=empty_msg(True)
            self.buf=self.cur.lines
            self.msg_buf.append(self.cur)
        self.buf[-1] += stuff
        return(len(stuff))
    def flush(self):
        if len(self.first.lines) > 0:
            self.ww.write_msg(self.first)
            self.first=empty_msg()
        if self.cur is not None:
            if self.cur.is_tick:
                self.count+=sum(map(len,self.buf))
            self.cur= None
            self.buf=None
        map(lambda TTT : self.ww.write_msg(TTT),self.msg_buf)
        self.msg_buf=list()
    def push(self):
        self.flush()
        self.last.lines.append('II: {:8d} items'.format(self.count))
        self.ww.write_msg(self.last)
        self.last=empty_msg()
        self.count=0


class MessageTicker:
    "Buffer to reduce amount of locking done on StdoutWordWrap"
    def __init__(self,ww):
        self.ww=ww
        self.count=0
        self.buf=None
    def write_msg(self,msg):
        if self.buf is not None:
            self.flush()
        self.ww.write_msg(msg)
        if msg.is_tick:
            tmp=sum(map(len,msg.lines))
            self.count+=tmp
            return(tmp)
        return(0)
    def write(self,stuff):
        if self.buf is not None:
            self.flush()
        self.ww.write_msg(PostBlock(False,[stuff]))
        return(0)
    def tick(self,stuff):
        if self.buf is None:
            self.buf=''
        self.buf += stuff
        if len(self.buf) > 100:
            self.flush()
        return(len(stuff))
    def flush(self):
        if self.buf is not None:
            self.count+=len(self.buf)
            self.ww.write_msg(PostBlock(True,self.buf))
            self.buf=None
    def push(self):
        self.flush()
        self.ww.write_msg(PostBlock(False,['II: {:8d} items'.format(self.count)]))
        self.count=0

class TickerWriter:
    "A suitable target for print to write ticks to message ticker"
    def __init__(self,t):
        self.t=t
    def write(self,stuff):
        self.t.tick(stuff)
