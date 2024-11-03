from bamRefine.functions import *
import os

dirN = os.path.dirname(os.path.realpath(__file__))
os.chdir(dirN)

lookup_l = 10
lookup_r = 10

### Double-stranded mode (default) -------------------------
snps_ds = "random_reads.snp"
snps_ds = parseSNPs(snps_ds, False)

def read_sam_file(filename):
    """Helper function to read SAM files"""
    stdout, _ = run_samtools('view', filename)
    for line in stdout.decode('utf-8').split('\n'):
        if line.strip():
            read_dict = parse_sam_line(line)
            if read_dict:
                yield read_dict

def test_flagReads_complex():
    inName = "complex_read.sam"
    for read_dict in read_sam_file(inName):
        mask, m_pos, m_side = flagReads(snps_ds, read_dict, lookup_l, lookup_r, read_dict)

        assert mask == "mask" 
        assert m_pos == [1,5,8,-9]
        assert m_side == [0,0,0,1]

def test_flagReads_complex_asymmetric():
    inName = "complex_read.sam"
    for read_dict in read_sam_file(inName):
        mask, m_pos, m_side = flagReads(snps_ds, read_dict, 10, 8, read_dict)

        assert mask == "mask" 
        assert m_pos == [1,5,8]
        assert m_side == [0,0,0]

def test_flagReads_short():
    inName = "short_read.sam"
    for read_dict in read_sam_file(inName):
        mask, m_pos, m_side = flagReads(snps_ds, read_dict, lookup_l, lookup_r, read_dict)
        print(m_pos)
        print(m_side)
        print(snps_ds)

        assert mask == "mask" 
        assert m_pos == [0,1,5,-2,-5]
        assert m_side == [0,0,0,1,1]

## ---------------------------------------------------

## single-stranded mode
snps_ss = "random_reads_ss.snp"
snps_ss = parseSNPs(snps_ss, singleStranded=True)

def test_flagReads_complex_single():
    inName = "complex_read.sam"
    for read_dict in read_sam_file(inName):
        mask, m_pos, m_side = flagReads(snps_ss, read_dict, lookup_l, lookup_r, read_dict)

        assert mask == "mask" 
        assert m_pos == [1,5,8]
        assert m_side == [0,0,0]

def test_flagReads_short_single():
    inName = "short_read.sam"
    for read_dict in read_sam_file(inName):
        mask, m_pos, m_side = flagReads(snps_ss, read_dict, lookup_l, lookup_r, read_dict)
        print(m_pos)
        print(m_side)
        print(snps_ds)

        assert mask == "mask" 
        assert m_pos == [0,1,5,-5,-9,-10]
        assert m_side == [0,0,0,1,1,1]