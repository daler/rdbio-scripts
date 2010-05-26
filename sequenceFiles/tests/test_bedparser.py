"""Test functions for bedparser.py"""

import bedparser
from cStringIO import StringIO
from nose import with_setup

def test_checkfields():
    """Test a single-line bedfile."""
    data = StringIO("""chrX\t10\t100\n""")
    a = bedparser.bedfile(data)
    a = list(a)
    i = a[0]
    assert i.chr=='chrX'
    assert i.start==10
    assert i.stop==100
    assert len(a) == 1

def test_9fields():
    'parse single.track.9.fields.bed'
    f = open('inputfiles/single.track.9.fields.bed')
    a = bedparser.bedfile(f)
    a = list(a)
    assert len(a) == 4
    check_9fields_first(a[0])
    check_9fields_last(a[-1])

    # check by passing a filename instead of open file
    f = 'inputfiles/single.track.9.fields.bed'
    a = bedparser.bedfile(f)
    a = list(a)
    assert len(a) == 4
    check_9fields_first(a[0])
    check_9fields_last(a[-1])

def test_singlefeature_filename():
    '''single-feature bed file with an extra newline on the end.'''
    f = 'inputfiles/single.track.unnamed.single.feature.bed'
    a = bedparser.bedfile(f)
    a = list(a)
    assert len(a) == 1
    check_singlefeature(a[0])

def check_9fields_first(i):
    'checks first feature in single.track.9.fields.bed'
    assert i.chr == 'chrX'
    assert i.start == 1 
    assert i.stop == 100
    assert i.name == 'feature1'
    assert i.value == 0.5
    assert i.strand == '+'
    assert i.thickStart == 1
    assert i.thickStop == 100
    assert i.itemRGB == '0,0,255'
    assert i.track.name is None
    assert i.track.description is None
    assert i.track.priority is None
    assert i.track.itemRgb is None
    assert i.track.visibility is None
    assert i.track.color is None
    assert i.track.offset is None
    assert i.track.url is None
    assert i.track.db is None
    assert i.track.group is None

def check_9fields_last(i):
    'checks last feature in single.track.9.fields.bed'
    assert i.chr == 'chrX'
    assert i.start == 5000 
    assert i.stop == 10000
    assert i.name == 'feature1'
    assert i.value == 90
    assert i.strand == '+'
    assert i.thickStart == 5000
    assert i.thickStop == 10000
    assert i.itemRGB == '0,0,255'
    assert i.track.name is None
    assert i.track.description is None
    assert i.track.priority is None
    assert i.track.itemRgb is None
    assert i.track.visibility is None
    assert i.track.color is None
    assert i.track.offset is None
    assert i.track.url is None
    assert i.track.db is None
    assert i.track.group is None

def check_singlefeature(i):
    'checks first feature in single.track.9.fields.bed'
    assert i.chr == 'chrX'
    assert i.start == 1 
    assert i.stop == 10
    assert i.name == None
    assert i.value == None
    assert i.strand == None
    assert i.thickStart == None
    assert i.thickStop == None
    assert i.itemRGB == None
    assert i.track.name is None
    assert i.track.description is None
    assert i.track.priority is None
    assert i.track.itemRgb is None
    assert i.track.visibility is None
    assert i.track.color is None
    assert i.track.offset is None
    assert i.track.url is None
    assert i.track.db is None
    assert i.track.group is None

def test_emptyfile():
    """Behavior upon getting an empty file"""

    # check with an empty, open file-like object.
    a = bedparser.bedfile(StringIO(''))
    a = list(a)
    assert len(a)==0
    
    # now check with a filename.
    a = bedparser.bedfile('inputfiles/empty.bed')
    a = list(a)
    assert len(a) == 0

def test_tracknames():
    """Bed file with a trackline definition"""
    a = bedparser.bedfile('inputfiles/multi.tracks.3.fields.bed')
    a = list(a)
    check_multitracks_first(a[0])
    check_multitracks_last(a[-1])

    # now check with an open file.
    a = bedparser.bedfile(open('inputfiles/multi.tracks.3.fields.bed'))
    a = list(a)
    check_multitracks_first(a[0])
    check_multitracks_last(a[-1])

def check_multitracks_first(i):
    'checks assertions on first record in multi.tracks.3.fields.bed'
    assert i.start == 1
    assert i.stop == 100
    assert i.chr == 'chrX'
    assert i.itemRGB is None
    assert i.value is None
    assert i.thickStart is None
    assert i.thickStop is None
    assert i.strand is None
    assert i.blockCount is None
    assert i.blockSizes is None
    
    # check track.
    assert i.track.name == '"track 1"'
    assert i.track.description == '"first track"'
    assert i.track.priority == 1
    assert i.track.itemRgb is None
    assert i.track.visibility is None
    assert i.track.color is None
    assert i.track.offset is None
    assert i.track.url is None
    assert i.track.db is None
    assert i.track.group is None

def check_multitracks_last(i):
    'checks assertions on last record in multi.tracks.3.fields.bed'
    assert i.start == 1
    assert i.stop == 9
    assert i.chr == 'chr3'
    assert i.itemRGB is None
    assert i.value is None
    assert i.thickStart is None
    assert i.thickStop is None
    assert i.strand is None
    assert i.blockCount is None
    assert i.blockSizes is None
    
    # check track.
    assert i.track.name == '"track 2"'
    assert i.track.description == '"second track"'
    assert i.track.priority is None
    assert i.track.itemRgb == 1
    assert i.track.visibility == 2
    assert i.track.color == '0,0,0'
    assert i.track.offset == 2
    assert i.track.url is None
    assert i.track.db == 'dm3'
    assert i.track.group is None

