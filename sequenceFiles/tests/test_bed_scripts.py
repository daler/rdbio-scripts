from sequenceFiles.bedSizeFilter import bedSizeFilter
import os,sys, tempfile
from cStringIO import StringIO
from nose import with_setup

tmp = tempfile.mktemp()

def test_bedSizeFilter_importing():
    infn = 'inputfiles/multi.tracks.3.fields.bed'
    outfn = tempfile.mktemp()
    infile = open(infn)
    outfile = open(outfn,'w')
    bedSizeFilter(infile,outfile,minlen=None,maxlen=None)
    outfile.close()
    infile.close()
    assert open('inputfiles/expected.bedSizeFilter.0.bed').read() == open(outfn).read()

def test_bedSizeFilter_cmdline():
    infn = 'inputfiles/multi.tracks.3.fields.bed'
    outfn = tempfile.mktemp()
    cmds = ' '.join(['python',
                     '../bedSizeFilter.py',
                     '-i', infn,
                     '-o', outfn,
                     ])
    os.system(cmds)
    assert open('inputfiles/expected.bedSizeFilter.0.bed').read() == open(outfn).read()

    # make sure minlen works
    infile = open(infn)
    outfile = open(outfn,'w')
    bedSizeFilter(infile,outfile,minlen=100,maxlen=None)
    outfile.close()
    assert open('inputfiles/expected.bedSizeFilter.1.bed').read() == open(outfn).read()

    # make sure minlen and maxlen work
    infile = open(infn)
    outfile = open(outfn,'w')
    bedSizeFilter(infile,outfile,minlen=10,maxlen=500)
    outfile.close()
    assert open('inputfiles/expected.bedSizeFilter.2.bed').read() == open(outfn).read()

def test_bedSizeFilter_pipe():
    infn = 'inputfiles/multi.tracks.3.fields.bed'
    outfn = tempfile.mktemp()
    # do the same test, but from the command line.
    # test stdin/stdout.
    cmds = ' '.join(['cat', infn,'|',
                     '../bedSizeFilter.py',
                     '--min','10',
                     '--max','500',
                     '>',outfn])
    os.system(cmds)
    assert open('inputfiles/expected.bedSizeFilter.2.bed').read() == open(outfn).read()




