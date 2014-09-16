#!/usr/bin/env python
'''A VCFv4.0 parser for Python.

The intent of this module is to mimic the ``csv`` module in the Python stdlib,
as opposed to more flexible serialization formats like JSON or YAML.  ``vcf``
will attempt to parse the content of each record based on the data types
specified in the meta-information lines --  specifically the ##INFO and
##FORMAT lines.  If these lines are missing or incomplete, it will check
against the reserved types mentioned in the spec.  Failing that, it will just
return strings.

There is currently one piece of interface: ``Reader``.  It takes a file-like
object and acts as a reader::

    >>> import vcf
    >>> vcf_reader = vcf.Reader(open('test/example.vcf', 'rb'))
    >>> for record in vcf_reader:
    ...     print record
    Record(CHROM=20, POS=14370, REF=G, ALT=['A'])
    Record(CHROM=20, POS=17330, REF=T, ALT=['A'])
    Record(CHROM=20, POS=1110696, REF=A, ALT=['G', 'T'])
    Record(CHROM=20, POS=1230237, REF=T, ALT=['.'])
    Record(CHROM=20, POS=1234567, REF=GTCT, ALT=['G', 'GTACT'])


This produces a great deal of information, but it is conveniently accessed.
The attributes of a Record are the 8 fixed fields from the VCF spec plus two
more.  That is:

    * ``Record.CHROM``
    * ``Record.POS``
    * ``Record.ID``
    * ``Record.REF``
    * ``Record.ALT``
    * ``Record.QUAL``
    * ``Record.FILTER``
    * ``Record.INFO``

plus three more attributes to handle genotype information:

    * ``Record.FORMAT``
    * ``Record.samples``
    * ``Record.genotypes``

``samples`` and ``genotypes``, not being the title of any column, is left lowercase.  The format
of the fixed fields is from the spec.  Comma-separated lists in the VCF are
converted to lists.  In particular, one-entry VCF lists are converted to
one-entry Python lists (see, e.g., ``Record.ALT``).  Semicolon-delimited lists
of key=value pairs are converted to Python dictionaries, with flags being given
a ``True`` value. Integers and floats are handled exactly as you'd expect::

    >>> vcf_reader = vcf.Reader(open('test/example.vcf', 'rb'))
    >>> record = vcf_reader.next()
    >>> print record.POS
    14370
    >>> print record.ALT
    ['A']
    >>> print record.INFO['AF']
    [0.5]

``record.FORMAT`` will be a string specifying the format of the genotype
fields.  In case the FORMAT column does not exist, ``record.FORMAT`` is
``None``.  Finally, ``record.samples`` is a list of dictionaries containing the
parsed sample column and ``record.genotypes`` is a dictionary of sample names
to genotype data::

    >>> record = vcf_reader.next()
    >>> for sample in record.samples:
    ...     print sample['GT']
    0|0
    0|1
    0/0
    >>> print record.genotypes['NA00001']['GT']
    0|0

Metadata regarding the VCF file itself can be investigated through the
following attributes:

    * ``VCFReader.metadata``
    * ``VCFReader.infos``
    * ``VCFReader.filters``
    * ``VCFReader.formats``
    * ``VCFReader.samples``

For example::

    >>> vcf_reader.metadata['fileDate']
    '20090805'
    >>> vcf_reader.samples
    ['NA00001', 'NA00002', 'NA00003']
    >>> vcf_reader.filters
    {'q10': Filter(id='q10', desc='Quality below 10'), 's50': Filter(id='s50', desc='Less than 50% of samples have data')}
    >>> vcf_reader.infos['AA'].desc
    'Ancestral Allele'

'''
import collections
import re
import csv

try:
    from ordereddict import OrderedDict
except ImportError:
    from collections import OrderedDict



# Metadata parsers/constants
RESERVED_INFO = {
    'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
    'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
    'END': 'Integer', 'H2': 'Flag', 'MQ': 'Float', 'MQ0': 'Integer',
    'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag', 'VALIDATED': 'Flag'
}

RESERVED_FORMAT = {
    'GT': 'String', 'DP': 'Integer', 'FT': 'String', 'GL': 'Float',
    'GQ': 'Float', 'HQ': 'Float'
}


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])


class _vcf_metadata_parser(object):
    '''Parse the metadat in the header of a VCF file.'''
    def __init__(self, aggressive=False):
        super(_vcf_metadata_parser, self).__init__()
        self.aggro = aggressive
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+)=(?P<val>.+)''')

    def read_info(self, info_string):
        '''Read a meta-information INFO line.'''
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError(
                "One of the INFO lines is malformed: {}".format(info_string))

        try:
            num = int(match.group('number'))
            if num < 0:
                num = None if self.aggro else '.'
        except ValueError:
            num = None if self.aggro else '.'

        info = _Info(match.group('id'), num,
                     match.group('type'), match.group('desc'))

        return (match.group('id'), info)

    def read_filter(self, filter_string):
        '''Read a meta-information FILTER line.'''
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: {}".format(
                    filter_string))

        filt = _Filter(match.group('id'), match.group('desc'))

        return (match.group('id'), filt)

    def read_format(self, format_string):
        '''Read a meta-information FORMAT line.'''
        match = self.format_pattern.match(format_string)
        if not match:
            raise SyntaxError(
                "One of the FORMAT lines is malformed: {}".format(
                    format_string))

        try:
            num = int(match.group('number'))
            if num < 0:
                num = None if self.aggro else '.'
        except ValueError:
            num = None if self.aggro else '.'

        form = _Format(match.group('id'), num,
                       match.group('type'), match.group('desc'))

        return (match.group('id'), form)

    def read_meta(self, meta_string):
        match = self.meta_pattern.match(meta_string)
        return match.group('key'), match.group('val')


# Reader class
class _meta_info(object):
    '''Decorator for a property stored in the header info.'''
    def __init__(self, func):
        self.func = func

    def __call__(self, fself):
        if getattr(fself, "_%s" % self.func.__name__) is None:
            fself._parse_metainfo()

        return self.func(fself)

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __doc__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__


class _Record(object):
    def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, genotypes):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        self.FILTER = FILTER
        self.INFO = INFO
        self.FORMAT = FORMAT
        self.genotypes = genotypes

    def __str__(self):
        return "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s)" % self.__dict__

    def add_format(self, fmt):
        self.FORMAT = self.FORMAT + ':' + fmt

    def add_filter(self, flt):
        if self.FILTER == '.':
            self.FILTER = ''
        else:
            self.FILTER = self.FILTER + ';'
        self.FILTER = self.FILTER + flt

    def add_info(self, info, value=True):
        self.INFO[info] = value

    @property
    def samples(self):
        """ return a list of samples, added for backwards compatibility """
        return self.genotypes.values()


class Reader(object):
    '''Read and parse a VCF v 4.0 file'''
    def __init__(self, fsock, aggressive=False):
        super(VCFReader, self).__init__()
        self.aggro = aggressive
        self._metadata = None
        self._infos = None
        self._filters = None
        self._formats = None
        self._samples = None
        self.reader = fsock
        self._header_lines = []
        if aggressive:
            self._mapper = self._none_map
        else:
            self._mapper = self._pass_map

    def __iter__(self):
        return self

    @property
    @_meta_info
    def metadata(self):
        '''Return the information from lines starting "##"'''
        return self._metadata

    @property
    @_meta_info
    def infos(self):
        '''Return the information from lines starting "##INFO"'''
        return self._infos

    @property
    @_meta_info
    def filters(self):
        '''Return the information from lines starting "##FILTER"'''
        return self._filters

    @property
    @_meta_info
    def formats(self):
        '''Return the information from lines starting "##FORMAT"'''
        return self._formats

    @property
    @_meta_info
    def samples(self):
        '''Return the names of the genotype fields.'''
        return self._samples

    def _parse_metainfo(self):
        '''Parse the information stored in the metainfo of the VCF.

        The end user shouldn't have to use this.  She can access the metainfo
        directly with ``self.metadata``.'''
        for attr in ('_metadata', '_infos', '_filters', '_formats'):
            setattr(self, attr, {})

        parser = _vcf_metadata_parser()

        line = self.reader.next()
        while line.startswith('##'):
            self._header_lines.append(line)
            line = line.strip()

            if line.startswith('##INFO'):
                key, val = parser.read_info(line)
                self._infos[key] = val

            elif line.startswith('##FILTER'):
                key, val = parser.read_filter(line)
                self._filters[key] = val

            elif line.startswith('##FORMAT'):
                key, val = parser.read_format(line)
                self._formats[key] = val

            else:
                key, val = parser.read_meta(line.strip())
                self._metadata[key] = val

            line = self.reader.next()

        fields = line.split()
        self._samples = fields[9:]

    def _none_map(self, func, iterable, bad='.'):
        '''``map``, but make bad values None.'''
        return [func(x) if x != bad else None
                for x in iterable]

    def _pass_map(self, func, iterable, bad='.'):
        '''``map``, but make bad values None.'''
        return [func(x) if x != bad else bad
                for x in iterable]

    def _parse_info(self, info_str):
        '''Parse the INFO field of a VCF entry into a dictionary of Python
        types.

        '''
        entries = info_str.split(';')
        retdict = {}
        for entry in entries:
            entry = entry.split('=')
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if entry_type == 'Integer':
                vals = entry[1].split(',')
                val = self._mapper(int, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = self._mapper(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type == 'String':
                val = entry[1]

            try:
                if self.infos[ID].num == 1:
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict

    def _parse_samples(self, samples, samp_fmt):
        '''Parse a sample entry according to the format specified in the FORMAT
        column.'''
        samp_data = OrderedDict()
        samp_fmt = samp_fmt.split(':')
        for name, sample in zip(self.samples, samples):
            sampdict = dict(zip(samp_fmt, sample.split(':')))
            for fmt in sampdict:
                vals = sampdict[fmt].split(',')
                try:
                    entry_type = self.formats[fmt].type
                except KeyError:
                    try:
                        entry_type = RESERVED_FORMAT[fmt]
                    except KeyError:
                        entry_type = 'String'

                if entry_type == 'Integer':
                    sampdict[fmt] = self._mapper(int, vals)
                elif entry_type == 'Float' or entry_type == 'Numeric':
                    sampdict[fmt] = self._mapper(float, vals)
                elif sampdict[fmt] == './.' and self.aggro:
                    sampdict[fmt] = None

            sampdict['name'] = name
            samp_data[name] = sampdict

        return samp_data

    def next(self):
        '''Return the next record in the file.'''
        if self._samples is None:
            self._parse_metainfo()
        row = self.reader.next().split()
        chrom = row[0]
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None if self.aggro else row[2]

        ref = row[3]
        alt = self._mapper(str, row[4].split(','))
        qual = 0 if row[5] == '.' or row[5] == "Infinity" else (float(row[5]) if '.' in row[5] else int(row[5]))
        filt = row[6].split(';') if ';' in row[6] else row[6]
        if filt == 'PASS' and self.aggro:
            filt = None
        info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None
            samples = None
        else:
            samples = self._parse_samples(row[9:], fmt)

        record = _Record(chrom, pos, ID, ref, alt, qual, filt, info, fmt,
                         samples)
        return record


class Writer(object):

    fixed_fields = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()

    def __init__(self, stream, template):
        self.writer = csv.writer(stream, delimiter="\t")
        self.template = template
        # TODO: shouldnt have to poke the parser to get the meta
        if template._samples is None:
            template._parse_metainfo()
        for line in template._header_lines:
            stream.write(line)
        self.write_header()

    def write_header(self):
        # TODO: write INFO, etc
        self.writer.writerow(self.fixed_fields + self.template.samples)

    def write_record(self, record):
        ffs = [record.CHROM, record.POS, record.ID, record.REF, self._format_alt(record.ALT),
            record.QUAL, record.FILTER, self._format_info(record.INFO), record.FORMAT]

        samples = [self._format_sample(record.FORMAT, sample)
            for sample in record.samples]
        self.writer.writerow(ffs + samples)

    def _format_alt(self, alt):
        return ','.join(alt)

    def _format_info(self, info):
        return ';'.join(["%s=%s" % (x, self._stringify(y)) for x, y in info.items()])

    def _format_sample(self, fmt, sample):
        if sample["GT"] == "./.":
            return "./."
        return ':'.join((str(self._stringify(sample[f])) for f in fmt.split(':')))

    def _stringify(self, x):
        if type(x) == type([]):
            return ','.join(map(str, x))
        return str(x)


def __update_readme():
    import sys
    file('README.rst', 'w').write(sys.modules[__name__].__doc__)


# backwards compatibility
VCFReader = Reader
VCFWriter = Writer

