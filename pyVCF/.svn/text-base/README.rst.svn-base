A VCFv4.0 parser for Python.

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

    * ``Reader.metadata``
    * ``Reader.infos``
    * ``Reader.filters``
    * ``Reader.formats``
    * ``Reader.samples``

For example::

    >>> vcf_reader.metadata['fileDate']
    '20090805'
    >>> vcf_reader.samples
    ['NA00001', 'NA00002', 'NA00003']
    >>> vcf_reader.filters
    {'q10': Filter(id='q10', desc='Quality below 10'), 's50': Filter(id='s50', desc='Less than 50% of samples have data')}
    >>> vcf_reader.infos['AA'].desc
    'Ancestral Allele'

Random access is supported for files with tabix indexes.  Simply call fetch for the
region you are interested in:

    >>> vcf_reader = vcf.Reader(filename='test/tb.vcf.gz')
    >>> for record in vcf_reader.fetch('20', 1110696-1, 1230237):
    ...     print record
    Record(CHROM=20, POS=1110696, REF=A, ALT=['G', 'T'])
    Record(CHROM=20, POS=1230237, REF=T, ALT=['.'])


