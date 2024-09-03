#!/usr/bin/env python

# Embedded file name: fill_blank.py
import sys, time
def process_file(input_fn, output_fn):
    try:
        input_file = open(input_fn, 'r')
    except:
        print 'IOError: [Errno 2] No such file or directory: %s' % input_fn

    try:
        output_file = open(output_fn, 'w')
    except:
        print 'IOError: [Errno 2] No permission for file directory: %s' % output_fn

    header = 0
    ZeroReadCount = 0
    RC = 0
    newline = ''
    for line in input_file:
        if line[0] == '@':
            if header == 1:
                output_file.write(newline)
                RC += 1
            header = 1
            newline = line
        elif line[0] == '\n' and header == 1:
            newline += 'NNNNN\n+\n$$$$$\n'
            output_file.write(newline)
            ZeroReadCount += 1
            header = 0
            RC += 1
            newline = ''
            continue
        else:
            newline += line

    if header == 1:
        output_file.write(newline)
        RC += 1
    input_file.close()
    output_file.close()
    return (RC, ZeroReadCount)


if __name__ == '__main__':
    start = time.time()
    try:
        RC, ZRC = zip(process_file(sys.argv[1], sys.argv[2]))
    except:
        print 'please insert 2 args, input file and output file name\n'
        print '--usage example\n'
        print '  python fill_blank.pyc input.fastq output.fastq '
        exit()

    print 'input file : %s' % sys.argv[1]
    print 'process time : ', time.time() - start, 's'
    print 'Total input reads : %s' % RC
    print 'Filling Blank reads : %s' % ZRC
    print 'output file : %s' % sys.argv[2]
