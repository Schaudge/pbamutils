pbamutils
============

The parallel version of some bam statistics utilities based on mmbam library


Introduction
------------
Partial bam information statistics functions include:
    read counts
    mpileup


Compile & Install
-----------------
Standard command to install the project:

```
git clone 
mkdir build && cd build
cmake .. && make
```

If some errors occur in compile the project, you can solve the problem by the following two ways:

1. replace the config.h file in dirctionary lib/mmbam, by the new created config.h filed follow the document https://yiq.gitlab.io/mmbam/ by (GNU) Autotools. And we support some code snippets for easy use in this project for reference:
```
cd lib/mmbam
autoreconf -i && ./configure
cp config.h src/mmbam
```

2. review and check the configurated value in config.h file


Extension
---------
Other utilities would by be extended easily following the source codes.

