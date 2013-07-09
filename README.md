bio-db-vcf
----------

DESCRIPTION
-----------

Easy parsing of Variant Call Format (VCF)

FEATURES/PROBLEMS
-----------------

* Works with format versions 4.1 and 4.0
* Support for earlier versions is not tested

SYNOPSIS
--------

```ruby
require 'bio'
require 'bio/db/vcf'

Bio::FlatFile.open(Bio::Db::Vcf, "example.vcf").each_entry do |e|
end
```

REQUIREMENTS
------------

The following gems are needed
* bio
* parslet

INSTALL
-------

FOR DEVELOPERS
--------------

After checking out the source, run:

```
$ rake newb
```

This task will install any missing dependencies, run the tests/specs,
and generate the RDoc.

LICENSE
-------

(The MIT License)

Copyright (c) 2013 Fedor Gusev

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
