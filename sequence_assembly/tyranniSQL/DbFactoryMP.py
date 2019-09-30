#!/usr/bin/env python2

from __future__ import print_function
import itertools, os
from multiprocessing import RLock, Pool, Process

"""

Name: PipelineDb.py

Author: Graham E Derryberry
Date: 11 September 2016

Description:
Tools for handling MP access to PipelineDb.

"""

global db_factory_singleton = None
db_pool = None;

def master_init(pdb):
	"Called to set connection arguments for the workers andd initialize the pool"
	pass

def worker_init(db_factory):
	"Passed to the pool to initialize the workers"
	pass


class DbFactory(Object):
	"""The factory object
	
	As I plan on making this some form of singleton, I probably should make this class
	a new style class that inherits from Object. 
	
	Upon initiation, the singleton instance either inherits or creates a PipelineDb 
	object, known as master.
	
	MP calls are worked such that when foo.pldb is None, DbFactory.getPLDB() is called.
	DbFactory then grabs the singleton, and checks to see if the current process has a 
	PipelineDb object. If not, master is cloned to create the current process instance.
	The instance for current process is then returned.
	
	"""
	
	def __new__(cls):
		"""Grab the singleton"""
		pass
	self._master
	self._local
	def getPLDB():
		"""Get the pipeline db for the current process"""
		if db_factory_singleton is None:
			db_factory_singleton=DbFactory()
		return db_factory_singleton._getPLDB()
	def _getPLDB(self):
		"""The local instance accessor"""
		pass