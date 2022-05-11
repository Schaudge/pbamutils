Mmbam README
============

Memory mapped parallel BAM file access API for high throughput sequence analysis informatics

Introduction
------------

High throughput, genome wide next generation sequencing (NGS) have revolutionized precision medicine. Increasingly, genomics-guided precision medicine is helping advanced cancer patients who have exhausted standard-of-care options. In these settings, the amount of data analyzed is small compared to large cohort studies, involving usually one to a few tumor samples and the paired normal sample from the same patient. However, the analysis turnaround is of critical importance. For example, it is the longstanding standard of care practice to initiate definitive therapy within 2 to 4 weeks after pediatric cancer patients receive biopsy / resection surgery. Furthermore, after the optimal treatment is identified, it still takes significant time to coordinate treatment access due to e.g. drug acquisition, compassionate care approval, clinical trial enrollment, or insurance authorization. It is therefore significant that the informatics analysis tasks, which have surpassed sequencing as the primary bottleneck, are to be as fast as current computer hardware can make possible.

We developed mmbam, which uses memory mapped I/O and the bam file index (BAI) for parallel data reading, and takes the scatter / gather programming paradigm to parallelize computation tasks over many different genomic regions. Since the majority of sequence analysis informatics analyses (e.g. quality control, various types of mutation calling) involve reading already indexed BAM files, mmbam has the potential to significantly shorten end-to-end analysis turnaround.

Getting started
---------------

We recommend that you download our pre-packaged distribution tarball at https://gitlab.com/yiq/mmbam/-/jobs/1655766656/artifacts/raw/build/mmbam-0.1.0.tar.gz since the repository does not contain autotools genenerated files (e.g. the configure script).

If you prefer to clone the repository, you need to install autoconf, automake, and libtools, as well as to run `autoreconf -i` in order to properly regenerate the build system.

Comprehensive documentation is available at https://yiq.gitlab.io/mmbam/

