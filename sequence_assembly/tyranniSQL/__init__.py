import argparse
import tyranniSQL.PipelineDb
import mysql.connector

def add_db_args(parser):
    parser.add_argument(
            "-u","--user",
            default="tyranni.user",
            type=str,
            help="""The username to use with the database [%(default)s]"""
        )
    parser.add_argument(
            "-p","--password",
            default="XXXXXX",
            type=str,
            help="""The password for the provided account"""
        )
    parser.add_argument(
            "-t","--host",
            default="eeb-globus.tulane.edu",
            type=str,
            help="""The MySQL database host [%(default)s]"""
        )
    parser.add_argument(
            "--port",
            default="3306",
            type=str,
            help="""The MySQL database port [%(default)s]"""
        )
    parser.add_argument(
            "--db",
            default="tyranni_pipeline",
            type=str,
            help="""The database to use [%(default)s]"""
        )
    parser.add_argument(
            "-n","--dry-run",
            action="store_true",
            help="""Only scan the file for likely actions, do not act on database""")
    parser.add_argument(
            "--no-commit",dest="rollback",
            action="store_true",
            help="""Insert, update and then roll back edits.  Will increase AUTO_INCREMENT silently.""")
    return parser

DEFAULT=argparse.Namespace(user="tyranni.user",
    password="XXXXXX",
    db="tyranni_pipeline",
    host="travailtu.mooo.com",
    port="43306",
    dry_run=True,
    rollback=True)

TEST=argparse.Namespace(user="tyranni.user",
    password="XXXXXX",
    db="test_pipeline",
    host="travailtu.mooo.com",
    port="43306",
    dry_run=True,
    rollback=True)

COMMIT=argparse.Namespace(user="tyranni.user",
    password="XXXXXX",
    db="tyranni_pipeline",
    host="travailtu.mooo.com",
    port="43306",
    dry_run=False,
    rollback=False)

def get_connection_from_args(args, delayed=False, ticker=None):
    "Makes Pipeline object based on the args from the provided parser"
    if delayed:
        cArgs={user:args.user,
            password:args.password,
            host:args.host,port:args.port,
            database:args.db}
        pdb=tyranniSQL.PipelineDb.PipelineDb(cxn_args=cArgs,
            rollback=args.rollback,
            dry_run=args.dry_run,ticker=ticker)
        return(pdb)
    conn= mysql.connector.connect(user=args.user,
        password=args.password,
        host=args.host,port=args.port,
        database=args.db)
    pdb=tyranniSQL.PipelineDb.PipelineDb(conn=conn,rollback=args.rollback,dry_run=args.dry_run,ticker=ticker)
    return(pdb)


# public symbols
#__all__ = [ "DEFAULT","TEST","COMMIT","get_connection_from_args","add_db_args" ]
