import argparse
import logging
import sys


parser = argparse.ArgumentParser()
parser.add_argument("bedpe")

logging.basicConfig(level=logging.INFO)

valid_chrom = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
    '11', '12', '13', '14', '15', '16', '17', '18', '19',
    '20', '21', '22', 'x', 'y', 'mt',
]
valid_strand = ['+', '-', '.']


class BedpeRow(object):

    def __init__(self, chrom1, start1, end1, chrom2, start2,
                 end2, name='', score='', strand1='.', strand2='.',
                 *extra):

        self.chrom1 = chrom1
        self.start1 = start1
        self.end1 = end1
        self.chrom2 = chrom2
        self.start2 = start2
        self.end2 = end2
        self.name = name
        self.score = score
        self.strand1 = strand1
        self.strand2 = strand2
        self.extra = extra


class InvalidChrom1(Exception): pass
class InvalidChrom2(Exception): pass
class InvalidStrand1(Exception): pass
class InvalidStrand2(Exception): pass
class InvalidStartEnd1(Exception): pass
class InvalidStartEnd2(Exception): pass


def validate_bedpe_row(row):
    if row.chrom1.lower() not in valid_chrom:
        raise InvalidChrom1()
    if row.chrom2.lower() not in valid_chrom:
        raise InvalidChrom2()
    if row.strand1 not in valid_strand:
        raise InvalidStrand1()
    if row.strand2 not in valid_strand:
        raise InvalidStrand2()
    if int(row.start1) + 1 > int(row.end1):
        raise InvalidStartEnd1()
    if int(row.start2) + 1 > int(row.end2):
        raise InvalidStartEnd2()


def fix_strand(strand):
    if strand == '1':
        return '+'
    elif strand == '-1':
        return '-'
    else:
        return strand


def fix_start_end(start, end):
    if int(start) + 1 > int(end):
        return end, start
    else:
        return start, end


def format_row(row):
    return '\t'.join([
        row.chrom1,
        row.start1,
        row.end1,
        row.chrom2,
        row.start2,
        row.end2,
        row.name,
        row.score,
        row.strand1,
        row.strand2,
    ] + list(row.extra))


def parse(path):
    rows = open(path).read().splitlines()
    outputs = []
    for i, line in enumerate(rows):
        sp = line.split("\t")
        row = BedpeRow(*sp)

        row.strand1 = fix_strand(row.strand1)
        row.strand2 = fix_strand(row.strand2)
        row.start1, row.end1 = fix_start_end(row.start1, row.end1)
        row.start1, row.end1 = fix_start_end(row.start1, row.end1)

        try:
            validate_bedpe_row(row)
        except InvalidChrom1:
            m = "Dropping row. Invalid chrom 1 on row {}: '{}'"
            m = m.format(i, row.chrom1)
            logging.info(m)
            continue
        except InvalidChrom2:
            m = "Dropping row. Invalid chrom 2 on row {}: '{}'"
            m = m.format(i, row.chrom2)
            logging.info(m)
            continue
        except InvalidStrand1:
            m = "Invalid strand 1 on row {}: '{}'"
            m = m.format(i, row.strand1)
            logging.error(m)
            sys.exit(1)
        except InvalidStrand2:
            m = "Invalid strand 2 on row {}: '{}'"
            m = m.format(i, row.strand2)
            logging.error(m)
            sys.exit(1)
        except InvalidStartEnd1:
            m = "Invalid start/end 1 on row {}: {} {}"
            m = m.format(i, row.start1, row.end1)
            logging.error(m)
            sys.exit(1)
        except InvalidStartEnd2:
            m = "Invalid start/end 2 on row {}: {} {}"
            m = m.format(i, row.start2, row.end2)
            logging.error(m)
            sys.exit(1)

        outputs.append(format_row(row))

    sys.stdout.write('\n'.join(outputs))

if __name__ == "__main__":
    args = parser.parse_args()
    parse(args.bedpe)
