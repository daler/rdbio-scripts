#!/usr/env/python
"""
This script converts a BED-format file to a sqlite3 database.

The database has the following schema:

CREATE TABLE features (chrom text,
                       start int,
                       stop int, 
                       name text,
                       value float,
                       strand text,
                       thickStart int,
                       thickStop int, 
                       itemRGB text,
                       blockCount int,
                       blockSizes text,
                       blockStarts text);
"""
import bedparser
import sqlite3
import optparse

op = optparse.OptionParser()
op.add_option('-i',dest='bedfile',help='Input bed file')
op.add_option('--db',dest='database',help='Output sqlite3 database')

def bed2db(bedfn, dbfn):
    conn = sqlite3.connect(dbfn)
    c = conn.cursor()
    c.execute('''
    CREATE TABLE features (chrom text,
                           start int,
                           stop int, 
                           name text,
                           value float,
                           strand text,
                           thickStart int,
                           thickStop int, 
                           itemRGB text,
                           blockCount int,
                           blockSizes text,
                           blockStarts text)
    ''')
    for feature in bedparser.bedfile(bedfn):
        items = [feature.chr, 
                 feature.start, 
                 feature.stop, 
                 feature.name, 
                 feature.value, 
                 feature.strand, 
                 feature.thickStart,
                 feature.thickStop, 
                 feature.itemRGB,
                 feature.blockCount,
                 feature.blockSizes, 
                 feature.blockStarts]
        c.execute('''
        INSERT INTO features VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
        ''',items)

    c.execute('CREATE INDEX idx_starts ON features(start)')
    c.execute('CREATE INDEX idx_stops ON features(stop)')
    c.execute('CREATE INDEX idx_chroms ON features(chrom)')
    c.execute('CREATE INDEX idx_values ON features(value)')
    conn.commit()
    return conn

if __name__ == "__main__":
    options,args = op.parse_args()
    bed2db(options.bedfn, options.dbfn)

