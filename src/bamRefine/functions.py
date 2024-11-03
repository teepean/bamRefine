import os
import pickle
import sys
import subprocess as sp
import json
from importlib.metadata import version

from .sam_parser import SamParser, SamReader, SamWriter

# Initialize SAM parser
sam_parser = SamParser()

def run_samtools(command, *args, input_data=None, capture_output=True):
    """Helper function to run samtools commands"""
    cmd = ['samtools'] + [command] + list(args)
    if input_data:
        process = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = process.communicate(input=input_data)
    else:
        process = sp.run(cmd, capture_output=capture_output)
        stdout = process.stdout
        stderr = process.stderr
    return stdout, stderr

def isDangerous(var):
    if ("C" in var) and ("G" in var):
        return (True, 'both')
    elif ("C" in var):
        return (True, '0')
    elif ("G" in var):
        return (True, '1')
    else:
        return (False, None)

def isDangerous_single(var):
    if ("C" in var):
        return (True, '0')
    else:
        return (False, None)

def generateTags(snpList, sideList):
    tag1 = ",".join([str(sideList.count(x)) for x in [0,1]])
    tag1 = ('ZC',tag1)

    sideSwitchI = "".join(str(x) for x in sideList).find('1')
    if sideSwitchI != -1:
        tag2_1 = ",".join([str(x) for x in snpList[0:sideSwitchI]])
        tag2_2 = ",".join([str(x) for x in snpList[sideSwitchI:]])
    else:
        tag2_1 = ",".join([str(x) for x in snpList])
        tag2_2 = ""
    tag2 = tag2_1 + ";" + tag2_2
    tag2 = ('ZP', tag2)

    return [tag1, tag2]

def modifyHeader(header_text, bamrefine_version, bamrefine_command_line):
    """Modify SAM header to add bamrefine program group"""
    header_lines = header_text.strip().split('\n')
    pg_entries = [line for line in header_lines if line.startswith('@PG')]
    
    # Find previous bamrefine entries
    prev_bamrefine = []
    for pg in pg_entries:
        fields = dict(field.split(':', 1) for field in pg.split('\t')[1:])
        if 'bamrefine' in fields.get('ID', ''):
            prev_bamrefine.append(fields['ID'])
    
    # Generate new ID
    this_bamrefine = ""
    if prev_bamrefine:
        if '.' in prev_bamrefine[-1]:
            last_bamrefine = prev_bamrefine[-1]
            this_bamrefine = last_bamrefine.split(".")[-1]
            this_bamrefine = str(int(this_bamrefine) + 1)
            this_bamrefine += "."
        else:
            this_bamrefine = ".1"
    ID = "bamrefine" + this_bamrefine
    
    # Add new PG line
    new_pg = f"@PG\tID:{ID}\tPN:bamrefine\tVN:{bamrefine_version}\t"
    if pg_entries:
        new_pg += f"PP:{pg_entries[-1].split('\t')[1].split(':')[1]}\t"
    new_pg += f"CL:{bamrefine_command_line}"
    
    # Insert new PG line before first non-header line
    for i, line in enumerate(header_lines):
        if not line.startswith('@'):
            header_lines.insert(i, new_pg)
            break
    else:
        header_lines.append(new_pg)
    
    return '\n'.join(header_lines)

def flagReads(snpLocDic, bamLine, look_l, look_r, read_dict):
    """
    Flag positions to be masked. Mapped read positions in <read_dict> that
    overlap with selected variants in the <snpLocDic> within <look_l> bases
    from the 5' end and <look_r> bases from the 3' end will be flagged for
    masking.
    """
    chrm = bamLine['ref_name']
    seq = bamLine['seq']
    snpList = [] # store positions to be masked
    sideList= [] # store mask sides of positions (5' or 3')

    ## make look_l = look_r if the latter not set:
    look_r = look_l if look_r is None else look_r
    look_list = [look_l, look_r]

    refpos = get_reference_positions(read_dict)
    refpos = [x + 1 if x is not None else x for x in refpos] ## convert to 1-based indices
    read_len = len(seq)

    ## Default inspectRange:
    inspectRange = [refpos[:look_l], refpos[read_len-1:read_len-look_r-1:-1]]

    ## Adjust if read is too short for the lookup values from either side:
    for side in range(2):
        sign = [1,-1][side]
        if read_len <= look_list[side]:
            inspectRange[side] = refpos[::sign] ## reverse the pos list if 3' end

    for side in range(2):
        for i, nt in enumerate(inspectRange[side]):
            shift = [0,1][side]
            sign = [1,-1][side]
            i = i + shift
            i = i * sign ### get indices for the 3' end
            key = chrm + " " + str(nt)
            try:
                snp = snpLocDic[side * -1][key] ## Always use the 5' side if ss-mode
            except KeyError:
                continue
            else:
                snpList.append(i)
                sideList.append(side)

    if len(snpList) > 0:
        return ('mask', snpList, sideList)
    else:
        return ('nomask', snpList, sideList)

def detectSNPfileFormat(fName):
    if fName.endswith('.snp'):
        chrI = 1
        posI = 3
        refI, altI = (4, 5)
    elif fName.endswith('.bed'):
        # Read first line and count fields using Python
        with open(fName) as f:
            first_line = f.readline()
            ncol = len(first_line.strip().split())

        chrI, posI, refI, altI = tuple(range(4))

        if ncol == 5:
           posI, refI, altI = tuple([x+1 for x in (posI, refI, altI)])

    return (chrI, posI, refI, altI)

def parseSNPs(fName, singleStranded):
    categorize_snps_func = isDangerous
    if singleStranded:
        categorize_snps_func = isDangerous_single

    snps = [{}, {}]

    chrI, posI, refI, altI = detectSNPfileFormat(fName)

    with open(fName) as snpF:
        for snp in snpF:
            snp = snp.strip().split()
            curC = snp[chrI]
            dangerous, side = categorize_snps_func([snp[refI], snp[altI]])
            if not dangerous:
                # ignoring otherwise
                continue
            key = curC + " " + snp[posI]
            if side != 'both':
                snps[int(side)][key] = [snp[x] for x in [chrI, posI, refI, altI]]
            else:
                snps[0][key] = [snp[x] for x in [chrI, posI, refI, altI]]
                snps[1][key] = [snp[x] for x in [chrI, posI, refI, altI]]

    if len(snps[1]) == 0:
        del(snps[1])

    return snps

def processBAM(inBAM_path, ouBAM_path, snps, contig, lookup, addTags=False):
    """Process BAM file to mask positions and add tags"""
    lookup = lookup.split(",")
    lookup_l = int(lookup[0])
    try:
        lookup_r = int(lookup[1])
    except IndexError:
        lookup_r = None

    stats = [0, 0]

    # Get header from input BAM
    stdout, _ = run_samtools('view', '-H', inBAM_path)
    header = stdout.decode('utf-8')

    # Process reads using optimized SAM parser
    with SamWriter(ouBAM_path, header=header) as writer:
        for read_dict in SamReader(inBAM_path, region=contig):
            mask, m_pos, m_side = flagReads(snps, read_dict, lookup_l, lookup_r, read_dict)
            
            if mask == 'mask':
                # Modify sequence and quality scores
                t = list(read_dict['seq'])
                q = list(read_dict['qual'])
                for p in m_pos:
                    t[p] = 'N'
                    q[p] = '!'
                read_dict['seq'] = "".join(t)
                read_dict['qual'] = "".join(q)
                
                # Update statistics
                st = [m_side.count(x) for x in [0,1]]
                for x in range(2):
                    stats[x] += st[x]
                    
                # Add tags if requested
                if addTags:
                    tags = generateTags(m_pos, m_side)
                    for tag_name, tag_value in tags:
                        read_dict['tags'][tag_name] = ('Z', tag_value)
                        
                writer.write(read_dict)
            elif mask == 'nomask':
                writer.write(read_dict)

    # Write statistics
    with open(contig+"_stats.txt", 'w') as statsF:
        statsF.write(str(stats[0]) + "\t" + str(stats[1]) + '\n')

def distributeSNPs(fName, chrms, singleStranded):
    categorize_snps_func = isDangerous
    if singleStranded:
        categorize_snps_func = isDangerous_single

    snps = {chrm: [{}, {}] for chrm in chrms}

    chrI, posI, refI, altI = detectSNPfileFormat(fName)

    with open(fName) as snpF:
        for snp in snpF:
            snp = snp.strip().split()
            curC = snp[chrI]
            dangerous, side = categorize_snps_func([snp[refI], snp[altI]])
            if not dangerous:
                # ignoring otherwise
                continue
            key = curC + " " + snp[posI]
            if side != 'both':
                snps[curC][int(side)][key] = [snp[x] for x in [chrI, posI, refI, altI]]
            else:
                snps[curC][0][key] = [snp[x] for x in [chrI, posI, refI, altI]]
                snps[curC][1][key] = [snp[x] for x in [chrI, posI, refI, altI]]

    ## IMPORTANT!!
    snps = {chrm: [snps[chrm][i] for i in range(len(snps[chrm])) if len(snps[chrm][i]) > 0] for chrm in snps if len(snps[chrm][0]) > 0}

    ss = ["ds", "ss"][int(singleStranded)]
    present_chrms = set()
    for chrm in snps:
        present_chrms.add(chrm)
        pickleName = os.path.basename(fName) + "_singleStranded_" + ss + "_contig_" + chrm +".brf"
        with open(pickleName, 'wb') as pfile:
            pickle.dump(snps[chrm], pfile)

    return present_chrms

def handleSNPs(fName, singleStranded, contig):
    ss = ["ds", "ss"][int(singleStranded)]
    pickleName = os.path.basename(fName) + "_singleStranded_" + ss + "_contig_" + contig +".brf"
    f_exists = os.path.isfile(pickleName)
    if f_exists:
        with open(pickleName, 'rb') as f:
            snps = pickle.load(f)
    return snps

def createBypassBED(inName, chrms, snps, singleStranded):
    s_chrms = distributeSNPs(snps, chrms, singleStranded)

    toBypass = set(chrms).difference(s_chrms)
    toFilter = s_chrms
    if toBypass:
        # Use samtools idxstats instead of pysam
        stdout, _ = run_samtools('idxstats', inName)
        bed = stdout.decode('utf-8').split('\n')
        bed = [x.split('\t') for x in bed if x]
        bed = [x[0:2] for x in bed if x[2] != '0']
        bed = [x for x in bed if x[0] in toBypass]
        bed = ["\t".join([x[0], "0", x[1]]) for x in bed]
        with open("bypass.bed", 'w') as bedF:
            for line in bed:
                bedF.write(line + '\n')
        bypass = True
    else:
        toFilter = chrms
        bypass = False

    return (list(toFilter), bypass)

def fetchChromosomes(inName):
    stdout, _ = run_samtools('idxstats', inName)
    idstats = stdout.decode('utf-8').split('\n')
    idstats = [x.split('\t') for x in idstats if x]
    chrms = [x[0] for x in idstats if x[2] != '0']
    return chrms

def main(args=None):
    inName = sys.argv[1]
    contig = sys.argv[2]
    lookup = sys.argv[3]
    snps = sys.argv[4]
    addTags = bool(int(sys.argv[5]))
    singleStranded = bool(int(sys.argv[6]))

    snps = handleSNPs(snps, singleStranded, contig)
    ouName = contig + ".bam"

    bamrefine_commandline = os.environ['BAMREFINE_CMDLINE']
    bamrefine_version = version("bamrefine")

    # Get header from input BAM
    stdout, _ = run_samtools('view', '-H', inName)
    header = stdout.decode('utf-8')
    
    # Modify header
    header = modifyHeader(header, bamrefine_version, bamrefine_commandline)

    # Process BAM file
    processBAM(inName, ouName, snps, contig, lookup, addTags)

    return 0