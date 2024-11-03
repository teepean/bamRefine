"""
Optimized SAM format parser module.
"""
import re
from typing import Dict, List, Tuple, Optional, Iterator
from dataclasses import dataclass
from collections import defaultdict

@dataclass
class CigarOp:
    """Represents a CIGAR operation"""
    length: int
    op: str

class CigarParser:
    """Fast CIGAR string parser with caching"""
    _cache = {}
    _cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    
    @classmethod
    def parse(cls, cigar: str) -> List[CigarOp]:
        """Parse CIGAR string into list of operations with caching"""
        if cigar in cls._cache:
            return cls._cache[cigar]
            
        ops = [CigarOp(int(length), op) 
               for length, op in cls._cigar_pattern.findall(cigar)]
        
        # Cache only if reasonable length to prevent memory issues
        if len(cigar) < 100:
            cls._cache[cigar] = ops
        return ops

    @classmethod
    def clear_cache(cls):
        """Clear the CIGAR operation cache"""
        cls._cache.clear()

class SamParser:
    """Fast SAM format parser"""
    
    def __init__(self):
        self.cigar_parser = CigarParser()
        
    def parse_line(self, line: str) -> Optional[Dict]:
        """Parse a single SAM line into a dictionary"""
        if not line or line.startswith('@'):
            return None
            
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11:
            return None
            
        # Parse tags efficiently
        tags = {}
        if len(fields) > 11:
            for tag in fields[11:]:
                try:
                    tag_name, tag_type, value = tag.split(':', 2)
                    tags[tag_name] = (tag_type, value)
                except ValueError:
                    continue
                    
        return {
            'name': fields[0],
            'flag': int(fields[1]),
            'ref_name': fields[2],
            'pos': int(fields[3]),
            'mapq': int(fields[4]),
            'cigar': fields[5],
            'rnext': fields[6],
            'pnext': int(fields[7]),
            'tlen': int(fields[8]),
            'seq': fields[9],
            'qual': fields[10],
            'tags': tags
        }

    def get_reference_positions(self, read_dict: Dict) -> List[Optional[int]]:
        """Get reference positions for a read based on CIGAR string"""
        pos = read_dict['pos'] - 1  # Convert to 0-based
        ref_positions = []
        seq_pos = 0
        
        cigar_ops = CigarParser.parse(read_dict['cigar'])
        
        for op in cigar_ops:
            if op.op in {'M', '=', 'X'}:  # Match or mismatch
                ref_positions.extend(range(pos, pos + op.length))
                pos += op.length
                seq_pos += op.length
            elif op.op in {'I', 'S'}:  # Insertion or soft clip
                ref_positions.extend([None] * op.length)
                seq_pos += op.length
            elif op.op in {'D', 'N'}:  # Deletion or skip
                pos += op.length
            elif op.op == 'H':  # Hard clip
                pass
            elif op.op == 'P':  # Padding
                pass
                
        return ref_positions

    def format_line(self, read_dict: Dict) -> str:
        """Convert a read dictionary back to SAM format"""
        fields = [
            read_dict['name'],
            str(read_dict['flag']),
            read_dict['ref_name'],
            str(read_dict['pos']),
            str(read_dict['mapq']),
            read_dict['cigar'],
            read_dict['rnext'],
            str(read_dict['pnext']),
            str(read_dict['tlen']),
            read_dict['seq'],
            read_dict['qual']
        ]
        
        for tag, (tag_type, value) in read_dict.get('tags', {}).items():
            fields.append(f"{tag}:{tag_type}:{value}")
        
        return '\t'.join(fields)

class SamReader:
    """Iterator class for reading SAM/BAM files using samtools"""
    
    def __init__(self, filename: str, region: Optional[str] = None):
        self.filename = filename
        self.region = region
        self.parser = SamParser()
        
    def __iter__(self) -> Iterator[Dict]:
        import subprocess as sp
        
        cmd = ['samtools', 'view', self.filename]
        if self.region:
            cmd.append(self.region)
            
        with sp.Popen(cmd, stdout=sp.PIPE, text=True) as proc:
            for line in proc.stdout:
                read_dict = self.parser.parse_line(line)
                if read_dict:
                    yield read_dict

class SamWriter:
    """Class for writing SAM/BAM files using samtools"""
    
    def __init__(self, filename: str, header: Optional[str] = None):
        self.filename = filename
        self.parser = SamParser()
        self.process = None
        self.header = header
        
    def __enter__(self):
        import subprocess as sp
        
        cmd = ['samtools', 'view', '-b', '-o', self.filename, '-']
        self.process = sp.Popen(cmd, stdin=sp.PIPE, text=True)
        
        if self.header:
            self.process.stdin.write(self.header)
        
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.process:
            self.process.stdin.close()
            self.process.wait()
            
    def write(self, read_dict: Dict):
        """Write a read dictionary as a SAM line"""
        if self.process and self.process.poll() is None:
            line = self.parser.format_line(read_dict)
            self.process.stdin.write(line + '\n')